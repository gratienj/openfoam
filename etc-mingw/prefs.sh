# ----------------------------------------------------------------------------
# Preferences for mingw cross-compilation

export FOAM_CONFIG_ETC="etc-mingw"

export WM_COMPILER=Mingw
export WM_MPLIB=msmpi-10.0

unset WM_COMPILE_CONTROL

# No zlib available:
# export WM_COMPILE_CONTROL="~libz"

# No mpi available:
# WM_MPLIB=none

# ----------------------------------------------------------------------------
