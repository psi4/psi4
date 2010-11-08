# Function to look for Apple vecLib on OS X systems
define(AC_CHECK_VECLIB,
[AC_PROVIDE([$0])
  echo "Entering ac_check_veclib"
  SAVE_LIBS=$LIBS
  echo "$target_vendor"
  if eval $CC
  LIBS=$SAVE_LIBS
  echo "Finished with ac_check_veclib"
]
)

