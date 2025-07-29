from run_builder import *

#~~~~~~~~~~~~~~~~~~~~~~~~ CUSTOM MATERIALS ~~~~~~~~~~~~~~~~~~~~~~~~#
# custom_material("ybco", 6.4, yttrium=13.35, barium=41.23, copper=28.62, oxygen=16.81)
# custom_material("Blade15Front", 8.8812266, copper=92.6, arsenic=2.20, silver=0.43, tin=0.68, antimony=0.43, lead=0.58, nickel=2.90, zinc=0.15)
# custom_material("Blade15Back",  8.8931427, copper=94.2, arsenic=1.98, silver=0.37, tin=0.61, antimony=0.19, lead=0.26, nickel=2.40)
#~~~~~~~~~~~~~~~~~~~~~~~~ CUSTOM MATERIALS ~~~~~~~~~~~~~~~~~~~~~~~~#

place("CarbonLayer", "cube", (0., 0., 0.), material= "Blade15Front", size_x = 10. * mm, size_y = 10. * mm, size_z = 90. * um)

place("source_marker", "sphere", (0., 0., -5. * cm), material="vacuum", radius = 2. * mm, alpha=0.5, red=0, green = 100)

place("detector", "sphere", (0., 0., 0.), material="vacuum", radius = 10 * cm, inner_radius = 9.5 * cm, alpha=0.01, red=255, green = 0)
make_sd("RBS", "detector", ["ekin", "theta", "phi"], "primary")
make_sd("PIXE", "detector", ["ekin", "theta", "phi"], "gamma")
make_beam_source("alpha", 4., (0, 0., -5. * cm), (0, 0, 1), sigma_r = 1* mm)

set_output_path("output_calsam1")
set_run_name("CalSam1")
config_run(8e9, 12) #Wenn Anwesend dann simulation sonst Interaktiv aber keine Simulation.
start_run("./output_runfiles/CalSam1_AlphaX.run") #Hier Pfad rein und der Code macht eine Datei mit dem Output.
