import os
import subprocess

# Create the virtual environment
subprocess.run(['python3.9', '-m', 'venv', 'myenv'])

subprocess.run('source myenv/bin/activate && pip install pandas', shell=True)

# Install packages inside the virtual environment
subprocess.run(['pip', 'install', 'pandas'])