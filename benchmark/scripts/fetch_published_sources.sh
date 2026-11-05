#!/usr/bin/env bash
set -euo pipefail

commit=316bf2a6bcd072336f4ed698776cad6b87f64db6
cache_dir="${1:-benchmark/cache/BayesPI-BAR}"
mkdir -p "$cache_dir"

curl -fL "https://raw.githubusercontent.com/junbaiw/BayesPI-BAR/${commit}/download/results/Mutations.tar.gz" -o "$cache_dir/Mutations.tar.gz"
curl -fL "https://raw.githubusercontent.com/junbaiw/BayesPI-BAR/${commit}/download/results/rankings.tar.gz" -o "$cache_dir/rankings.tar.gz"

printf '%s  %s\n' \
  9108a14efab606898e0c9c8d94fac1f842b717d79ac9294c919c83ac1805158d "$cache_dir/Mutations.tar.gz" \
  cf8d614dcd6e98be10af737a09820c6e18def0dd0be253d362d58755a226ae0a "$cache_dir/rankings.tar.gz" | sha256sum -c -

