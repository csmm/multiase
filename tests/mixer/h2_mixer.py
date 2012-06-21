from ase import Atoms
from gpaw import GPAW

from csmmcalc.mixer.mixer import Mixer
from csmmcalc.lammps.reaxff import ReaxFF
from csmmcalc.utils import get_datafile

d = 0.74
a = 6.0

atoms = Atoms("H2",
                positions = [(0, 0, 0),
                (0, 0, d)],
                cell = (10*a, 10*a, 10*a),
                tags=[Mixer.tag(0, 0.5) |
                    Mixer.tag(1, 0.5),
                    Mixer.tag(0, 0.5) |
                    Mixer.tag(1, 0.5)])

atoms.center()
calc_1 = GPAW(nbands=2, txt="h2_1.txt")
calc_2 = GPAW(nbands=2, txt="h2_2.txt")
calc_3 = ReaxFF(specorder = ("C", "H", "O", "N", "S"),
       ff_file_path=get_datafile("ffield.reax"))
#calc_3.keep_tmp_files = True
mixer = Mixer([calc_1, calc_3])
atoms.set_calculator(mixer)
print(atoms.get_forces())

