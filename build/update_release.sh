set -e
#@TODO: might want to change version here automatically
python increase_version.py
python set_test_mode_false.py
bash make_zip.sh
transcriptic upload-release release.zip "DelveTX2"
bash publish_public_protocols.sh
