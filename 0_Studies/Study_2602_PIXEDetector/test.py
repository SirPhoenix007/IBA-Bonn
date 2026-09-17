import argparse

parser = argparse.ArgumentParser(description='PIXE Spectrum Plotting')

parser.add_argument('-sm','--simple_mode', choices=['None','single', 'dual'], default='None', help='Plotting mode for simple plots: single or dual axis')
parser.add_argument('-dm','--detailed_mode', choices=['None','1','2'], default='None', help='Plotting mode for detailed plots: single axis + expected lines (1,2)')

parser.add_argument('-p','--path', help='Directory of data to insert into routine.')
parser.add_argument('-s','--save', help='Save the generated plot as .png and .pdf files', action='store_true')

if __name__ == "__main__":
    args = parser.parse_args()
    
    print(args.save)