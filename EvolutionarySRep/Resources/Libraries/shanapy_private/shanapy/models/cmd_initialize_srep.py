import vtk, qt, ctk, slicer
import sys
img_file_path = sys.argv[1]
output_file_path = sys.argv[2]
slicer.modules.skeletalrepresentationinitializer.logic().CLIInitialize(img_file_path, output_file_path, 5, 9)
print("Finished initialization ", output_file_path)
exit()