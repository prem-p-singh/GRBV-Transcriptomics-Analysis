#!/usr/bin/env python
"""00. Build per-sample weather-window aggregates.

Reads config sampling dates + CIMIS daily files, produces one CSV with
per-sample tavg_mean / tmax_max / gdd_sum / precip_sum / rh_mean for a
configurable window ending on each sample's sampling date.
"""
import argparse, pandas as pd, os, sys, yaml
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _weather_utils import per_sample_weather

ap = argparse.ArgumentParser()
ap.add_argument("--meta", required=True)
ap.add_argument("--weather_2021", required=True); ap.add_argument("--weather_2022", required=True)
ap.add_argument("--config", required=True); ap.add_argument("--out", required=True)
args = ap.parse_args()

os.makedirs(os.path.dirname(args.out), exist_ok=True)
meta = pd.read_csv(args.meta, sep="\t")
df = per_sample_weather(meta,
                        config_path=args.config,
                        weather_paths={2021: args.weather_2021, 2022: args.weather_2022})
df.to_csv(args.out)
print(f"wrote {len(df)} samples × {df.shape[1]} weather features")
