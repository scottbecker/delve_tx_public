set -e
rm release.zip
cd  ../src
zip -r ../build/release.zip *   --exclude=test* --exclude=manual_runs*
cd ../build
#must come last because we are overwritting manifest with different parameters
zip -j release.zip requirements.txt /tmp/manifest.json
rm /tmp/manifest.json
