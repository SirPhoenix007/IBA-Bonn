from run_builder import *

#~~~~~~~~~~~~~~~~~~~~~~~~ CUSTOM MATERIALS ~~~~~~~~~~~~~~~~~~~~~~~~#
# density in g/cm³
# custom_material("ybco", 6.4, yttrium=13.35, barium=41.23, copper=28.62, oxygen=16.81)
# custom_material("Blade15Front", 8.881, copper=92.6, arsenic=2.20, silver=0.43, tin=0.68, antimony=0.43, lead=0.58, nickel=2.90, zinc=0.15)
# custom_material("Blade15Back",  8.893, copper=94.2, arsenic=1.98, silver=0.37, tin=0.61, antimony=0.19, lead=0.26, nickel=2.40)
custom_material("Columbian_Coffee", 0.350, 
                carbon=81.6733, 
                oxygen=11.1815,
                nitrogen=4.3753,
                potassium=2.4231,
                phosphorus=0.573,
                magnesium=0.1371,
                calcium=0.0305,
                sodium=0.0079)
#~~~~~~~~~~~~~~~~~~~~~~~~ CUSTOM MATERIALS ~~~~~~~~~~~~~~~~~~~~~~~~#

place("target", "cube", (0., 0., 0.), material= "Columbian_Coffee", size_x = 10. * mm, size_y = 10. * mm, size_z = 100. * um)

place("source_marker", "sphere", (0., 0., -5. * cm), material="vacuum", radius = 2. * mm, alpha=0.5, red=0, green = 100)

place("detector", "sphere", (0., 0., 0.), material="vacuum", radius = 10 * cm, inner_radius = 9.5 * cm, alpha=0.1, red=255, green = 0)
make_sd("RBS", "detector", ["ekin", "theta", "phi"], "primary")
make_sd("PIXE", "detector", ["ekin", "theta", "phi"], "gamma")

beam_type = "proton"
energy = 2. #in MeV
make_beam_source(beam_type, energy, (0, 0., -5. * cm), (0, 0, 1), sigma_r = 1* mm)

set_output_path(f"Col_Coffee_{beam_type}_{str(energy).replace('.','-')}")
set_run_name(f"Col_Coffee_{beam_type}_{str(energy).replace('.','-')}_")
config_run(1e8, 12) #Wenn Anwesend dann simulation sonst Interaktiv aber keine Simulation.
start_run(f"./output_runfiles/Col_Coffee_{beam_type}_{str(energy).replace('.','-')}.run", force_file_write=True) #Hier Pfad rein und der Code macht eine Datei mit dem Output.
