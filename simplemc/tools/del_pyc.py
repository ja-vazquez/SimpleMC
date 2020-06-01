
# Short script to remove .pyc files (which dont want on github)

import os

root =  "/Users/josevazquezgonzalez/Desktop/Codigos_SimpleMC/SimpleMC"

for path, subdirs, files in os.walk(root):
    for name in files:
        if '.pyc' in name:
            print (os.path.join(path, name))
            os.remove(os.path.join(path, name))