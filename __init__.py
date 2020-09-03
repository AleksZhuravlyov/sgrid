import sys
import os

current_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(current_path))

from sgrid_bind import Sgrid
from sgrid_bind import save_files_collection_to_file