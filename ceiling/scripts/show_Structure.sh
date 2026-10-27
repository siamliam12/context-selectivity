python3 -c "
import os
def draw_tree(path, prefix=''):
    if not os.path.isdir(path): return
    files = sorted(os.listdir(path))
    for i, file in enumerate(files):
        if file.startswith('.'): continue # Skip hidden files
        is_last = (i == len(files) - 1)
        print(f\"{prefix}{'└── ' if is_last else '├── '}{file}\")
        if os.path.isdir(os.path.join(path, file)):
            draw_tree(os.path.join(path, file), prefix + ('    ' if is_last else '│   '))
draw_tree('.')
"
