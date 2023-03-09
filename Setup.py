import subprocess
import sys

def create_venv():
    subprocess.check_call([sys.executable, '-m', 'venv', 'venv'])
    subprocess.check_call(['./venv/bin/pip', 'install', '-r', 'requirements.txt'])

# if __name__ == '__main__':
#     create_venv()