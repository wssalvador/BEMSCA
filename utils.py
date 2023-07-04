import os, shutil

def clear_folders():
    folders = ['simul', 'workflows', 'results', 'studies']
    for folder in folders:
        for filename in os.listdir(folder):
            file_path = os.path.join(folder, filename)
            if os.path.isfile(file_path):
                os.unlink(file_path)