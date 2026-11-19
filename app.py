from __future__ import annotations

import csv
import io
import re
import zipfile
from datetime import date, datetime, timedelta
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import pandas as pd
import streamlit as st
from astroquery.jplhorizons import Horizons

HCT = {"lon": 78.96, "lat": 32.78, "elevation": 4500}
MONTHS = {
    "Jan": "01", "Feb": "02", "Mar": "03", "Apr": "04",
    "May": "05", "Jun": "06", "Jul": "07", "Aug": "08",
    "Sep": "09", "Oct": "10", "Nov": "11", "Dec": "12",
}

st.set_page_config(page_title="HCT Comet Trackfile Creator", page_icon="🔭", layout="wide")
st.title("🔭 Comet Trackfile Creator for HCT, Hanle")
st.caption("Generate HCT telescope trackfiles using ephemerides retrieved directly from NASA/JPL Horizons.")
st.info("Enter the observing date in UT. The observing month is a separate input used in the HCT trackfile header: make sure it matches the month in the date.")

c1, c2 = st.columns(2)
with c1:
    month = st.selectbox("Select observing month", list(MONTHS.keys()), index=datetime.utcnow().month - 1)
    observing_date = st.date_input("Date of observation (UT)", value=date.today())
    object_type = st.selectbox("Object type", ["Comet", "Asteroid"])
    number_objects = st.selectbox("Number of objects", list(range(1, 11)), index=0)
with c2:
    make_plot = st.checkbox("Generate visibility plot (includes the Moon)", value=True)
    st.markdown("**Object names must be recognised by JPL Horizons.**")
    st.caption("Examples: 3I/ATLAS, 12P/Pons-Brooks, C/2023 A3, 29P, 99942 Apophis")

st.markdown("---")
object_inputs = []
for i in range(number_objects):
    st.markdown(f"**Object {i + 1}**")
    a, b = st.columns([3, 1])
    with a:
        obj_name = st.text_input("Object name / designation", key=f"obj_name_{i}", placeholder="e.g. 12P/Pons-Brooks", label_visibility="collapsed")
    with b:
        disintegrated = st.checkbox("Disintegrated", key=f"dis_{i}", help="Use the alternate JPL Horizons query without id_type.")
    object_inputs.append((obj_name.strip(), disintegrated))

submitted = st.button("Generate trackfiles", type="primary", use_container_width=True)


def clean_name(name: str) -> str:
    return re.sub(r"[ /]", "", name)


def retrieve_track_eph(name: str, disintegrated: bool, start: datetime, stop: datetime, idtype: str):
    kwargs = dict(id=name, location=HCT,
                  epochs={"start": str(start.date()), "stop": str(stop.date()), "step": "5760"})
    if not disintegrated:
        kwargs["id_type"] = idtype
    return Horizons(**kwargs).ephemerides(
        quantities="4,5", closest_apparition=True,
        skip_daylight=False, airmass_lessthan=3.0
    )


def first_continuous_block(datetime_strings):
    line_start = 0
    for j in range(len(datetime_strings) - 1):
        try:
            if int(datetime_strings[j + 1][12:14]) - int(datetime_strings[j][12:14]) > 1:
                line_start = j + 1
                break
        except (ValueError, IndexError, TypeError):
            continue
    line_end = len(datetime_strings)
    for j in range(line_start, len(datetime_strings) - 1):
        try:
            if int(datetime_strings[j + 1][12:14]) - int(datetime_strings[j][12:14]) > 1:
                line_end = j + 1
                break
        except (ValueError, IndexError, TypeError):
            continue
    return line_start, line_end


def make_trackfile(name: str, eph, start: datetime, month_label: str):
    """Preserve the original GUI's table handling and HCT rate/azimuth conversions."""
    clean = clean_name(name)
    dt_strings = list(eph["datetime_str"])
    line_start, line_end = first_continuous_block(dt_strings)
    if line_end <= line_start:
        raise ValueError("JPL Horizons returned no usable continuous ephemeris block.")

    # Equivalent to the original temporary text file -> pandas column rewrite.
    rows = list(zip(
        dt_strings[line_start:line_end], eph["AZ"][line_start:line_end],
        eph["EL"][line_start:line_end], eph["AZ_rate"][line_start:line_end],
        eph["EL_rate"][line_start:line_end]
    ))
    raw = io.StringIO()
    writer = csv.writer(raw, delimiter=" ")
    raw.write("Date_(UT) HR:MN:SC Azi Elev dAZ*cosE d(ELV)/dt\n")
    writer.writerows(rows)
    # The desktop GUI removed CSV quote characters before splitting the date and time fields.
    raw = io.StringIO(raw.getvalue().replace('"', '').replace("'", ""))
    # Match the original split-on-space parsing and column interpretation.
    data = pd.read_csv(raw, sep=" ")
    data.columns = ["Date", "Time", "Azi", "Elev", "Azi_cosE", "Elev_rate"]

    zeroes = np.zeros(len(data.Date), dtype=np.int64)
    num_1 = np.arange(0, len(data.Date), 1) * 30
    az_rate = [x / np.cos(np.deg2rad(y)) for x, y in zip(data.Azi_cosE, data.Elev)]
    az_rate_final = [x / (60 * 60 * 60) for x in az_rate]
    el_rate_final = [x / (60 * 60 * 60) for x in data.Elev_rate]

    date_1 = re.sub("-", "/", str(data.Date.iloc[0])).replace(month_label, MONTHS[month_label])
    lines = [f"{date_1} {str(data.Time.iloc[0])[:-4]}", str(len(data.Date))]
    for values in zip(
        num_1, np.round(data.Azi, 5), np.round(az_rate_final, 5),
        zeroes, zeroes, np.round(data.Elev, 5), np.round(el_rate_final, 5),
        zeroes, zeroes, zeroes, zeroes, zeroes, zeroes
    ):
        # First write the intermediate .trk row, then apply the original awk azimuth correction.
        az = float(values[1])
        if az > 325:
            az -= 360
        lines.append(
            f"{int(values[0])}\t{az:8.5f}\t{float(values[2]):8.5f}\t0\t0\t"
            f"{float(values[5]):8.5f}\t{float(values[6]):8.5f}\t0\t0\t0\t0\t0\t0"
        )
    return clean, "\n".join(lines) + "\n"


def visibility_data(name: str, disintegrated: bool, start: datetime, idtype: str):
    stop = start + timedelta(days=1)
    kwargs = dict(id=name, location=HCT,
                  epochs={"start": str(start.date()), "stop": str(stop.date()), "step": "10m"})
    if not disintegrated and name != "301":
        kwargs["id_type"] = idtype
    eph = Horizons(**kwargs).ephemerides(
        quantities="4,5", closest_apparition=True, skip_daylight=False
    )
    ut_list = list(eph["datetime_str"])
    el_list = list(eph["EL"])
    line_num = 0
    for j in range(len(ut_list) - 1):
        try:
            if int(ut_list[j + 1][12:14]) - int(ut_list[j][12:14]) > 1:
                line_num = j + 1
        except (ValueError, IndexError, TypeError):
            continue
    times_block = ut_list[line_num:]
    elev_block = el_list[line_num:]
    if not times_block:
        raise ValueError("No visibility ephemeris points were returned.")
    if times_block[-1].endswith("00:00"):
        times_block[-1] = times_block[-1][:-5] + "23:59"
    dt_list = []
    aligned_elev = []
    for tstr, elev in zip(times_block, elev_block):
        parts = tstr.split()
        if len(parts) < 2:
            continue
        try:
            dt = datetime.strptime(parts[0] + " " + parts[1], "%Y-%m-%d %H:%M")
        except ValueError:
            try:
                dt = datetime.combine(start.date(), datetime.strptime(parts[1], "%H:%M").time())
            except ValueError:
                continue
        dt_list.append(dt)
        aligned_elev.append(float(elev))
    if not dt_list:
        raise ValueError("Could not parse the visibility times returned by JPL Horizons.")
    return {"name": name, "dt": dt_list, "elev": aligned_elev}


def create_visibility_figure(entries, start_date):
    fig, ax = plt.subplots(layout="tight", figsize=(12, 8))
    for entry in entries:
        style = "r--" if entry["name"] == "301" else "-"
        label = "Moon" if entry["name"] == "301" else entry["name"]
        ax.plot(mdates.date2num(entry["dt"]), entry["elev"], style, label=label)
    night = mdates.date2num(datetime.combine(start_date, datetime.strptime("13:30", "%H:%M").time()))
    plot_start = mdates.date2num(datetime.combine(start_date, datetime.strptime("08:00", "%H:%M").time()))
    plot_end = mdates.date2num(datetime.combine(start_date + timedelta(days=1), datetime.strptime("00:00", "%H:%M").time()))
    ax.axhline(20, xmin=0, xmax=1, ls="--", color="k")
    ax.axvline(night, ymin=0, ymax=1, ls="--", color="k")
    ax.text(night, 75, "Night", rotation=90, fontsize=20)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="both", width=1)
    ax.tick_params(axis="both", which="major", length=10, color="k")
    ax.tick_params(axis="both", which="minor", length=6, color="k")
    ax.set_xlabel(start_date.strftime("%Y-%m-%d") + " UT", fontsize=22)
    ax.set_ylabel("Elevation (deg)", fontsize=22)
    ax.tick_params(axis="x", which="both", direction="in", bottom=True, top=True, labelsize=20)
    ax.tick_params(axis="y", which="both", direction="in", left=True, right=True, labelsize=20)
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_color("black")
    ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=120))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.set_xlim(plot_start, plot_end)
    ax.set_ylim(7, 90)
    ax.legend(fontsize=20, loc="upper right", frameon=True)
    ax.set_title("Visibility curve", fontsize=24)
    return fig


if submitted:
    errors = []
    valid_objects = []
    if observing_date.month != list(MONTHS.values()).index(MONTHS[month]) + 1:
        st.error("The selected observing month does not match the month in the observation date. Please correct it before continuing.")
    elif any(not name for name, _ in object_inputs):
        st.error("Please enter a name/designation for every object.")
    else:
        valid_objects = object_inputs

    if valid_objects:
        start_dt = datetime.combine(observing_date, datetime.min.time())
        stop_dt = start_dt + timedelta(days=2)
        idtype = "designation" if object_type == "Comet" else "name"
        date_tag = start_dt.strftime("%Y_%m_%d")
        date_tag_compact = start_dt.strftime("%Y%m%d")
        files_out = {}
        visibility_entries = []
        progress = st.progress(0, text="Starting JPL Horizons queries…")
        status = st.status("Generating HCT trackfiles", expanded=True)

        for i, (name, disintegrated) in enumerate(valid_objects):
            status.write(f"Processing {name} (disintegrated: {'Yes' if disintegrated else 'No'})")
            try:
                eph = retrieve_track_eph(name, disintegrated, start_dt, stop_dt, idtype)
                clean, content = make_trackfile(name, eph, start_dt, month)
                filename = f"{clean}_{date_tag}.trk"
                files_out[filename] = content.encode("utf-8")
                status.write(f"Created {filename}")
            except Exception as exc:
                msg = f"Trackfile generation failed for {name}: {type(exc).__name__}: {exc}"
                errors.append(msg)
                files_out[f"{clean_name(name)}_{date_tag_compact}_jplerror.log"] = (msg + "\n").encode("utf-8")
                status.write("⚠️ " + msg)
            progress.progress((i + 1) / (len(valid_objects) + (len(valid_objects) + 1 if make_plot else 0)))

        if make_plot:
            # Original workflow also queried the Moon (JPL Horizons ID 301).
            for i, (name, disintegrated) in enumerate(valid_objects + [("301", True)]):
                label = "Moon" if name == "301" else name
                status.write(f"Retrieving visibility ephemeris for {label}")
                try:
                    visibility_entries.append(visibility_data(name, disintegrated, start_dt, idtype))
                except Exception as exc:
                    msg = f"Visibility query failed for {label}: {type(exc).__name__}: {exc}"
                    errors.append(msg)
                    status.write("⚠️ " + msg)
                progress.progress((len(valid_objects) + i + 1) / (len(valid_objects) + len(valid_objects) + 1))

        fig = None
        if visibility_entries:
            try:
                fig = create_visibility_figure(visibility_entries, observing_date)
                st.subheader("Visibility plot")
                st.pyplot(fig, use_container_width=True)
                plot_buffer = io.BytesIO()
                fig.savefig(plot_buffer, format="png", bbox_inches="tight", dpi=120)
                files_out[f"Visibility_{start_dt.strftime('%Y-%m-%d')}_UT.png"] = plot_buffer.getvalue()
                plt.close(fig)
            except Exception as exc:
                errors.append(f"Visibility plotting failed: {type(exc).__name__}: {exc}")
        elif make_plot:
            st.warning("No visibility curves could be generated.")

        status.update(label="Processing finished", state="complete", expanded=False)
        progress.progress(1.0, text="Finished")

        trackfiles = [name for name in files_out if name.endswith(".trk")]
        if trackfiles:
            st.success(f"Generated {len(trackfiles)} trackfile(s).")
            st.subheader("Download trackfiles")
            cols = st.columns(min(3, len(trackfiles)))
            for i, filename in enumerate(trackfiles):
                with cols[i % len(cols)]:
                    st.download_button(
                        f"⬇️ {filename}", data=files_out[filename],
                        file_name=filename, mime="text/plain", key=f"dl_{filename}"
                    )
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zf:
                for filename, data in files_out.items():
                    zf.writestr(filename, data)
            st.download_button(
                "📦 Download all outputs (ZIP)", data=zip_buffer.getvalue(),
                file_name=f"HCT_Trackfiles_{date_tag}.zip",
                mime="application/zip", type="primary", use_container_width=True
            )
        else:
            st.error("No trackfiles were generated. Check the errors below and the JPL Horizons object designation.")

        if errors:
            with st.expander(f"Processing messages / errors ({len(errors)})", expanded=True):
                for error in errors:
                    st.warning(error)

st.markdown("---")
st.caption("HCT: longitude 78.96° E, latitude 32.78° N, elevation 4500 m. Ephemerides are queried live from NASA/JPL Horizons.")
st.caption("Software workflow adapted from the original HCT Trackfile Created by Dr. Aravind Krishnakumar | University of Liège | ")
