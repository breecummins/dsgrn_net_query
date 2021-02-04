"""
Installs the package requirements and Git repositories needed for the DSGRN design interface
"""
import os, subprocess


def rmdir(dir):
    # change permissions of all files in dir to ensure automatic deletion
    for root, dirs, files in os.walk(dir):
        for d in dirs:
            os.chmod(os.path.join(root, d), 0o777)
        for f in files:
            os.chmod(os.path.join(root, f), 0o777)

    # completely delete dir and all its contents
    if dir[-1] == os.sep: dir = dir[:-1]
    files = os.listdir(dir)
    for file in files:
        if file == '.' or file == '..': continue
        path = dir + os.sep + file
        if os.path.isdir(path):
            rmdir(path)
        else:
            os.unlink(path)
    os.rmdir(dir)

# Clone the necessary repositories
os.makedirs("packages", exist_ok=True)
os.chdir("packages")
for dir in ['dsgrn_utilities']:
    if os.path.isdir(dir):
        rmdir(dir)

os.system('git clone https://github.com/breecummins/dsgrn_utilities')
os.chdir('dsgrn_utilities/')
subprocess.call(['bash install.sh'], shell=True)
os.chdir('tests/')
subprocess.call(['pytest'], shell=True)
os.chdir("../..")
