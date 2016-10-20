cp ../src/test/{test_bacteria_utils.py,test_custom_connection.py,test_custom_protocol.py,test_protocols.py,test_utils.py} ../../dtx_public/src/test
rsync -av --delete ../src/test/data ../../dtx_public/src/test --exclude *.pyc
rsync -av --delete ../build ../../dtx_public/ --exclude *.pyc
rsync -av --delete ../src/lib ../../dtx_public/src/ --exclude *.pyc
rsync -av --delete ../src/transcriptic_tools ../../dtx_public/src/ --exclude *.pyc
python sync_public_protocols.py
python build_public_manifest.py
(cd ../../dtx_public/ && git add -A && git commit -m "update public protocols" && git push)
