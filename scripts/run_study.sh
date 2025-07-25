#!/bin/bash
#SBATCH --job-name=mesa_param_study
#SBATCH --output=../logs/mesa_%j.log
#SBATCH --time=7-00:00:00

# Paths
MESA_DIR="/home/zx/Documents/mesa-r22.11.1"
BASE_DIR="/home/zx/Documents/Zh_Research/AGN_Star_8/parameter_study"
TEMPLATE="${BASE_DIR}/MESA_workdir"
RESULTS="${BASE_DIR}/results"
LOGS="${BASE_DIR}/logs"
SUMMARY="${LOGS}/summary.log"

# Parameter lists
densities=(1e-18 1e-17 5e-16 1e-16 5e-15 1e-15)
sound_speeds=(1e6)

# Prepare directories & summary
mkdir -p "$RESULTS" "$LOGS"
: > "$SUMMARY"   # truncate summary log

# Loop over parameters (sequential)
for rho in "${densities[@]}"; do
  for cs in "${sound_speeds[@]}"; do
    run_name="rho_${rho//./_}_cs_${cs//./_}"
    run_dir="${RESULTS}/${run_name}"
    
    echo "=== Starting ${run_name} ===" >> "$SUMMARY"
#     mkdir -p "$run_dir"
#     cp -r "${TEMPLATE}/." "$run_dir"
    
    # Remove extra output dirs if they exist
    for d in grid_png LOGS photos; do
      [ -d "${run_dir}/${d}" ] && rm -rf "${run_dir}/${d}"
    done

    # Update inlist_project but preserve trailing "!..." notes
    sed -i -E "s#^(x_ctrl\(1\) =)[^!]*(\!.*)?#\1 ${rho} \2#"   "$run_dir/inlist_project"
    sed -i -E "s#^(x_ctrl\(2\) =)[^!]*(\!.*)?#\1 ${cs}cd \2#" "$run_dir/inlist_project"

    # Run MESA
    pushd "$run_dir" >/dev/null
      ./clean >/dev/null 2>&1
      ./mk > $run_dir/run.log 2>&1
      ./rn > rn.log   2>&1
    popd >/dev/null

#     # Check for completion
#     if grep -q "termination code:" "$run_dir/run.log"; then
#       echo "${run_name}: Success" >> "$SUMMARY"
#     else
#       echo "${run_name}: Failed"  >> "$SUMMARY"
#     fi

  done
done

echo "All runs done. Summary in ${SUMMARY}."


# #!/bin/bash
# #SBATCH --job-name=mesa_param_study
# #SBATCH --output=../logs/mesa_%A_%a.log
# #SBATCH --time=7-00:00:00

# export MESA_DIR="/home/zx/Documents/mesa-r22.11.1"

# # Define parameter ranges
# # densities=(1e-18 1e-17 1e-16 1e-15)
# densities=(1e-18 1e-17 5e-16 1e-16 5e-15 1e-15)
# sound_speeds=(1e6)

# # Base directory setup
# BASE_TEMPLATE = /home/zx/Documents/Zh_Research/AGN_Star_8/parameter_study
# #BASE_DIR=$(dirname $(pwd))  # Move up one level from scripts directory
# MESA_TEMPLATE="${BASE_DIR}/MESA_workdir"
# RESULTS_DIR="${BASE_DIR}/results"
# LOGS_DIR="${BASE_DIR}/logs"

# # Create required directories
# mkdir -p "${RESULTS_DIR}" "${LOGS_DIR}"

# # Main execution loop
# for density in "${densities[@]}"; do
#   for sound_speed in "${sound_speeds[@]}"; do
#     # Create unique run identifier
#     RUN_NAME="rho_${density//./_}_cs_${sound_speed//./_}"
#     RUN_DIR="${RESULTS_DIR}/${RUN_NAME}"
    
#     echo "Starting run: ${RUN_NAME}"
#     echo "============================================"
    
#     # Create run directory and copy template
#     mkdir -p "${RUN_DIR}"
#     cp -r "${MESA_TEMPLATE}/"* "${RUN_DIR}"

#     # Modify inlist parameters
#     sed -i "/x_ctrl(1) =/s/= .*/= ${density}/" "${RUN_DIR}/inlist_project"
#     sed -i "/x_ctrl(2) =/s/= .*/= ${sound_speed}cd/" "${RUN_DIR}/inlist_project"
    
#     # Execute MESA commands
#     (
#       cd "${RUN_DIR}"
#       ./clean >/dev/null 2>&1
#       ./mk > build.log 2>&1 && ./rn > run.log 2>&1
      
#       # Validate completion
#       if grep -q "termination code:" run.log; then
#         echo "${RUN_NAME}: Success" >> "${LOGS_DIR}/summary.log"
#       else
#         echo "${RUN_NAME}: Failed" >> "${LOGS_DIR}/summary.log"
#       fi
#     ) &
    
#   done
# done

# # Wait for all background processes
# wait

# echo "Parameter study complete. Results summary:"
# cat "${LOGS_DIR}/summary.log"