import subprocess
import os
import shutil

def run_julia_scripts(scripts):
    """
    Execute a list of Julia scripts in order.
    """
    for script in scripts:
        print(f"Running Julia script: {script}")
        try:
            subprocess.run(["julia", script], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error occurred while running {script}: {e}")
            exit(1)

def move_data_files(source_folder, target_folder):
    """
    Move all data files from the source folder to the target folder.
    """
    print(f"Moving data files from {source_folder} to {target_folder}...")
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)  # Create the target folder if it doesn't exist

    for file_name in os.listdir(source_folder):
        source_path = os.path.join(source_folder, file_name)
        target_path = os.path.join(target_folder, file_name)

        if os.path.isfile(source_path):
            print(f"Moving {source_path} to {target_path}")
            shutil.move(source_path, target_path)


def run_python_scripts(scripts):
    """
    Execute a list of Python scripts in order.
    """
    for script in scripts:
        print(f"Running Python script: {script}")
        try:
            subprocess.run(["python3", script], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error occurred while running {script}: {e}")
            exit(1)



def main():
    # Define paths for data and scripts
    data_source_folder = os.path.join("data")
    data_target_folder = os.path.join("data/singlerun")

    julia_scripts = [
        os.path.join("code", "run_model.jl"),
        os.path.join("code", "experiments.jl"),
        os.path.join("code", "run_model_ntimes.jl"),
        os.path.join("code", "run_model_shock.jl"),
    ]

    python_scripts = [
        os.path.join("results", "main_plots.py"),
        os.path.join("results", "SF_plotting.py"),
    ]

    # Step 1: Run the first Julia script
    print("Step 1: Running the first Julia script...")
    run_julia_scripts([julia_scripts[0]])

    # Step 2: Move data files to the singlerun folder
    print("Step 2: Moving data files to singlerun folder...")
    move_data_files(data_source_folder, data_target_folder)

    # Step 3: Run the second Julia script
    print("Step 3: Running the second Julia script...")
    run_julia_scripts([julia_scripts[1]])
    move_data_files(data_source_folder, os.path.join("data/OFAT"))

    # Step 4: Run the third Julia script
    print("Step 4: Running the third Julia script...")
    run_julia_scripts([julia_scripts[2]])
    move_data_files(data_source_folder, os.path.join("data/multirun"))

    run_julia_scripts([julia_scripts[3]])
    move_data_files(data_source_folder, os.path.join("data/priceshocks"))

    # Step 5: Run Python scripts
    print("Step 5: Running Python scripts...")
    run_python_scripts(python_scripts)

    print("Pipeline executed successfully. Results are ready!")

    print("Pipeline executed successfully. Results are ready!")

if __name__ == "__main__":
    main()
