#class to get file name

import os
import sys

class FileName:

    def __init__(self):
        pass


    def get_filename(self, ending):
        filenames = []
        for i in os.listdir(os.getcwd()):
            if i.endswith(ending):
                filenames.append(i)

            else:
                continue
        filenames.sort()
        self.listfn=filenames
