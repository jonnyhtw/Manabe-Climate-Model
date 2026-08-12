#!/bin/bash
cd "$HOME/test_MCM/run_dir"

echo "=== Linking atmos_r30.nc from restart/initial to atmos_r30_in.nc ==="
ln -sf "$HOME/test_MCM/restart/initial/atmos_r30.nc" atmos_r30_in.nc
ln -sf "$HOME/test_MCM/restart/initial/atmos_r30.nc" atmos_r30.nc

ls -l atmos_r30_in.nc

echo "=== Running readr_clean ==="
echo "  71" > fort.5
./readr_clean < fort.5

echo "=== Results ==="
echo "--- irt.data ---"
[ -f irt.data ] && cat irt.data || echo "irt.data missing"

echo "--- iyear.data ---"
[ -f iyear.data ] && cat iyear.data || echo "iyear.data missing"
