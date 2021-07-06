# find files with a given tag in a given path

#import stuff
import os

def find_files(path, tag):
  file_paths = []
  for root, dirs, files in os.walk(path):
    for file in files:
      if tag in file:
        path = os.path.join(root,file)
        file_paths.append(path)
  return file_paths 



def find_dirs(path, tag):
  dir_paths = []
  for root, dirs, files in os.walk(path):
    for dir in dirs:
      if tag in dir:
        path = os.path.join(root, dir)
        dir_paths.append(path)
  return dir_paths 
