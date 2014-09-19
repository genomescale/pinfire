#!/usr/bin/python2.7

import os
import sys
import subprocess

loganalyser_binary = sys.argv[1]
log_paths = sys.argv[2:]

burnin = 36110000

for log_file_path in log_paths:
	folder_path, file_name = os.path.split(log_file_path)
	base = file_name[file_name.find("_"):file_name.find(".")]

	ess_file_name = "ess" + base + ".csv"
	ess_file_path = os.path.join(folder_path, ess_file_name)

	loganalyser_cmd = [loganalyser_binary, "-burnin", str(burnin), log_file_path, ess_file_path]

	try:
		loganalyser_result = subprocess.check_output(loganalyser_cmd)
		print(0)
		print(loganalyser_result)
	except CalledProcessError as loganalyser_error:
		print(loganalyser_error.return_code)
		print(loganalyser_error.output)
