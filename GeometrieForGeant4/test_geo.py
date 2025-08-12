from run_builder import *

custom_material("ybco", 6.4, yttrium=13.35, barium=41.23, copper=28.62, oxygen=16.81)
place("target", "cube", (0., 0., 0.), material= "ybco", size_x = 10. * mm, size_y = 10. * mm, size_z = 5. * um)

place("source_marker", "sphere", (0., 0., -5. * cm), material="vacuum", radius = 2. * mm, alpha=0.5, red=0, green = 100)

place("detector", "sphere", (0., 0., 0.), material="vacuum", radius = 10 * cm, inner_radius = 9.5 * cm, alpha=0.01, red=255, green = 0)
make_sd("RBS", "detector", ["ekin", "theta", "phi"], "primary")
make_sd("PIXE", "detector", ["ekin", "theta", "phi"], "gamma")
make_beam_source("alpha", 28., (0, 0., -5. * cm), (0, 0, 1), sigma_r = 0.4* mm)

set_output_path("test_geo_output")
set_run_name("custom_name")
#config_run(1e6, 8) #Wenn Anwesend dann simulation sonst Interaktiv aber keine Simulation.
start_run() #Hier Pfad rein und der Code macht eine Datei mit dem Output.
