If running into compiler issues on Mac try:

* clean anaconda un-install:
https://docs.anaconda.com/anaconda/install/uninstall/

* sort the xcode developer issues described here
https://discourse.mc-stan.org/t/compiling-error-after-upgrading-to-macos-big-sur/19309/3
by:
```
sudo rm -rf /Library/Developer/CommandLineTools
xcode-select --install
```

* re-install Anconda:
https://docs.anaconda.com/anaconda/install/mac-os/
