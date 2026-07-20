#!/bin/bash

# Test file to verify the various argument passing strategies.
# 1. Call the small driver_cli_cfg_fusion.py script without any options 
# 2. Call the small driver_cli_cfg_fusion.py script with a CLI options redirecting the output
# 3. Call the small driver_cli_cfg_fusion.py script with another config file, overloading the WR algo
# 4. Call the small driver_cli_cfg_fusion.py script with another config file, overloading the WR algo and redirecting the output

# Execute the python command.
python driver_cli_cfg_fusion.py
python driver_cli_cfg_fusion.py --experiment_name results/with_cli_options
python driver_cli_cfg_fusion.py --cfg ../data/cfg_test.ini
python driver_cli_cfg_fusion.py --cfg ../data/cfg_test.ini --experiment_name results/with_cli_options_and_config_file