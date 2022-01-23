class FileName:
	
	def __init__(self):
		self

	def get_filename(self, ending):
		
		import os, sys
		filenames = []
		for i in os.listdir(os.getcwd()):
			if i.endswith(ending):
				filenames.append(i)
			else:
				continue
				
		filenames.sort()	
		self.listfn=filenames
