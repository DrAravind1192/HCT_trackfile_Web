# HCT Trackfile Creator — Web App

A Streamlit web interface for the HCT trackfile workflow. It queries NASA/JPL Horizons, creates `.trk` files, optionally plots visibility in the page (including the Moon), and lets users download individual trackfiles or a ZIP of all outputs.

## Run locally

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt
streamlit run app.py
```

## Deploy free with Streamlit Community Cloud

1. Create a new GitHub repository (or a new branch in your existing repository).
2. Put `app.py` and `requirements.txt` in the repository root.
3. Go to https://share.streamlit.io/ and sign in with GitHub.
4. Choose **Create app**, select the repository and branch, and set the main file to `app.py`.
5. Deploy. Streamlit will install `requirements.txt` and provide a shareable `*.streamlit.app` URL.

## Notes

- No Tkinter, desktop display, shell `awk`, or local persistent output folder is required. Downloads are generated in memory for each run.
- The observation month remains a separate input, as in the desktop GUI, but the app checks that it matches the selected date before processing.
- JPL Horizons queries require outbound internet access. Queries can take time; free hosting resources and availability are not guaranteed.
- The application retains the original Horizons query settings, HCT site coordinates, 5760-minute trackfile ephemeris step, 10-minute visibility step, rate conversions, azimuth correction (>325 degrees subtract 360), and visibility plot limits/markers.
- Review each generated file before using it at the telescope.
