"""This adds .. to the system path so that this file can find the jobs folder."""
from pathlib import Path
import sys

DOTDOT = Path(__file__).resolve().parents[1]
sys.path.append(str(DOTDOT))
