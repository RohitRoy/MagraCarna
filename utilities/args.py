import os
import sys

from ..magracarna.motifs import MotifAssigner

class ArgParser(object):

	@classmethod
	def get_sources(cls, folder, ids, *extensions):
		filepaths = list()
		folderlist = os.listdir(folder)
		for extension in extensions:
			extension_length = len(extension)
			for filename in folderlist:
				if filename.endswith(extension):
					if not (ids and filename[:-extension_length] not in ids):
						filepaths.append(os.path.join(folder, filename))
		return sorted(filepaths)

	@classmethod
	def basename(cls, filepath, delimiters):
		basename = os.path.basename(filepath)
		for delimiter in delimiters:
			basename = basename.split(delimiter)
			basename = basename[:-1] if len(basename) > 1 else basename[:1]
			basename = delimiter.join(basename)
		return basename

	@classmethod
	def map_sources(cls, sourcesA, sourcesB):
		sourcesA = {cls.basename(file, "._"): file for file in sourcesA}
		sourcesB = {cls.basename(file, "._"): file for file in sourcesB}
		common = sorted(set(sourcesA).intersection(sourcesB))
		sources = list()
		for sourceid in common:
			sources.append((sourcesA[sourceid], sourcesB[sourceid]))
		return sources

	@classmethod
	def get_filepaths(cls, paths):
		filepaths = list()
		for path in paths:
			if os.path.isdir(path):
				folderfiles = list()
				for filename in os.listdir(path):
					filepath = os.path.join(path, filename)
					if filename.split("_")[0].isdigit():
						order_by = int(filename.split("_")[0])
						folderfiles.append((order_by, filepath))
					else:
						folderfiles.append((-1, filepath))
				filepaths += [path for _, path in sorted(folderfiles)]
			else:
				filepaths.append(path)
		return filepaths

	@classmethod
	def get_assigner(cls, motifpaths, mode="all"):
		motifsfiles = cls.get_filepaths(motifpaths)
		assigner = MotifAssigner.extract(*motifsfiles, mode=mode)
		return assigner

	@classmethod
	def read_idslist(cls, idslist):
		ids = None
		if idslist and os.path.exists(idslist):
			with open(idslist, 'r') as idsfile:
				ids = [line.strip() for line in idsfile.readlines()]
		return ids

