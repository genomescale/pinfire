#!/usr/bin/python2.7

import os

rerun = set(['08r_32s_16l_01k_02i_00p', '09r_16s_16l_01k_02i_00p', '09r_16s_16l_04k_02i_03p', '09r_16s_16l_04k_04i_03p', '09r_16s_16l_04k_02i_01p', '09r_16s_16l_02k_04i_00p', '09r_16s_16l_04k_04i_01p', '09r_16s_04i_08k_16l_07p', '09r_16s_04i_08k_16l_03p', '09r_16s_04i_08k_08l_07p', '09r_16s_08l_02k_04i_01p', '09r_16s_04i_08k_16l_00p', '09r_16s_08l_04k_04i_03p', '09r_16s_04l_04k_04i_02p', '09r_16s_04i_08k_08l_00p', '09r_16s_08l_04k_04i_01p', '09r_16s_04i_04k_16l_03p', '09r_16s_08l_04k_02i_02p', '09r_16s_04l_04k_04i_00p', '09r_16s_08l_04k_02i_01p', '09r_16s_16l_02k_02i_00p', '09r_16s_04i_02k_16l_00p', '09r_16s_16l_04k_02i_00p', '09r_16s_08l_01k_04i_00p', '09r_16s_16l_02k_04i_01p', '09r_16s_16l_04k_04i_02p', '09r_16s_16l_04k_02i_02p', '09r_16s_04i_08k_16l_02p', '09r_16s_04i_08k_16l_04p', '09r_16s_04i_08k_08l_06p', '09r_16s_08l_04k_04i_00p', '09r_16s_08l_02k_04i_00p', '09r_16s_04i_08k_16l_01p', '09r_16s_04i_08k_08l_05p', '09r_16s_04i_08k_16l_05p', '09r_16s_04i_08k_08l_01p', '09r_16s_08l_04k_04i_02p', '09r_16s_04i_04k_16l_02p', '09r_16s_16l_01k_04i_00p', '09r_16s_04l_04k_04i_01p', '09r_16s_16l_04k_04i_00p', '09r_16s_08l_02k_02i_00p', '09r_16s_04i_08k_16l_06p', '09r_16s_16l_02k_02i_01p'])

old_folder = "starbeast"
new_folder = "rerun"

for r in rerun:
	r_filename = r + ".xml"
	old_xml_path = os.path.join(old_folder, r_filename)
	new_xml_path = os.path.join(new_folder, r_filename)
	os.rename(old_xml_path, new_xml_path)