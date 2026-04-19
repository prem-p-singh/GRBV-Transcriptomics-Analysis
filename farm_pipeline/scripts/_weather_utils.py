"""Shared helper: per-sample CIMIS weather aggregates around sampling dates."""
import pandas as pd, numpy as np, yaml

def load_cimis(path):
    w = pd.read_csv(path)
    w.columns = [c.strip() for c in w.columns]
    w["Date"] = pd.to_datetime(w["Date"], format="%m/%d/%Y", errors="coerce")
    w = w.dropna(subset=["Date"])
    # Standard column names
    ren = {"Max Air Temp (C)": "tmax",
           "Min Air Temp (C)": "tmin",
           "Avg Air Temp (C)": "tavg",
           "Precip (mm)": "precip",
           "ETo (mm)": "eto",
           "Sol Rad (W/sq.m)": "solrad",
           "Avg Rel Hum (%)": "rh"}
    for k, v in ren.items():
        if k in w.columns:
            w[v] = pd.to_numeric(w[k], errors="coerce")
    return w[["Date"] + [v for v in ren.values() if v in w.columns]]

def window_agg(w, end_date, days_before, base_gdd=10):
    """Mean tavg, max tmax, GDD sum, precip sum in (end_date - days_before, end_date]."""
    end   = pd.Timestamp(end_date)
    start = end - pd.Timedelta(days=days_before)
    s = w[(w["Date"] > start) & (w["Date"] <= end)]
    if len(s) == 0:
        return {}
    gdd = np.clip(s.get("tavg", pd.Series(dtype=float)) - base_gdd, 0, None).sum()
    out = {
        "tavg_mean":   float(s["tavg"].mean())   if "tavg"   in s else np.nan,
        "tmax_max":    float(s["tmax"].max())    if "tmax"   in s else np.nan,
        "tmin_min":    float(s["tmin"].min())    if "tmin"   in s else np.nan,
        "gdd_sum":     float(gdd),
        "precip_sum":  float(s["precip"].sum())  if "precip" in s else np.nan,
        "rh_mean":     float(s["rh"].mean())     if "rh"     in s else np.nan,
    }
    return out

def per_sample_weather(meta, config_path, weather_paths):
    """Return DataFrame indexed by Sample with weather-window features.

    weather_paths = {2021: <csv>, 2022: <csv>}
    """
    cfg = yaml.safe_load(open(config_path))
    dates = {k: pd.to_datetime(v) for k, v in cfg["sampling_dates"].items()}
    win   = cfg["weather_window"]["days_before"]
    base  = cfg["weather_window"]["base_temp_gdd"]
    weather = {int(yr): load_cimis(p) for yr, p in weather_paths.items()}

    rows = []
    for _, r in meta.iterrows():
        key = f"{r['Sampling_Season']}|{r['Maturity_stage']}"
        if key not in dates:
            rows.append({"Sample": r["Sample"]})
            continue
        yr = int(r["Sampling_Season"])
        agg = window_agg(weather[yr], dates[key], win, base_gdd=base)
        agg["Sample"] = r["Sample"]
        rows.append(agg)
    return pd.DataFrame(rows).set_index("Sample")
