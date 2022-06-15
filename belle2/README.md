mkdir -p include/framework/logging
find /scratchfs/bes/myliu/cepc/CEPCSW_bak/CEPCSW/belle2/ -name "*.h" -exec ln -s {}  include/framework/logging \;
