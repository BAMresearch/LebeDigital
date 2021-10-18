import os
from pathlib import Path
import git

print('hello')
baseDir = Path(__file__).resolve().parents[0]
if os.path.exists(os.path.join(baseDir,'bamChowlk')):
    print('chowlk is already existed')
else:
    git.Git(baseDir).clone('https://github.com/firmao/bamChowlk')