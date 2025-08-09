python
import sys
import os
import glob

# Auto-detect GCC's pretty-printer path
matches = glob.glob("/usr/share/gcc-*/python")
if matches:
    pp_path = matches[0]
    sys.path.insert(0, pp_path)
    try:
        from libstdcxx.v6 import printers
        printers.register_libstdcxx_printers(None)
        print(f"[GDB] Loaded libstdc++ pretty-printers from {pp_path}")
    except Exception as e:
        print(f"[GDB] Failed to load pretty-printers: {e}")
else:
    print("[GDB] No libstdc++ pretty-printer path found")
end

set print pretty on
set print array on
set print elements 0
set print object on
set print static-members on
set print demangle on
set history save on
set pagination off
