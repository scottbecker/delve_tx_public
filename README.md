# delve_tx_public

This is a set of transcriptic package protocols.

## Protocol Overview

Bacteria Freezing: Take bacteria from a well, amplify it, freeze (-80C) 10 tubes w/ 115ul each of it
Bacteria Pelleting: Take bacteria from a well, amplify it, possibly induce to high copy number, and pellet
DNA Resuspension: Add TE to a tube with DNA, mix, measure concentration, and optionally split the tube (-80C) into a 4C storage working stock
Create Water Container: Handy test protocol to mimick reagent containers in test mode by making a container with every well filled to max volume with water



## Running protocols from the command line (requires unix)


1. update protocols/the_protocol/sample_input.json with the correct container id's you want to use
2. update ~/.transcriptic to have your login credentials

{
  "api_root": "https://secure.transcriptic.com",
  "token": "your token",
  "email": "your email",
  "organization_id": "your org"
}

3. Add the following to src/auth.json

{
  "X_User_Email":"your email",
  "X_User_Token":"your user token",
  "org_name": "your org name"
}

4. Add the following to src/test_mode_auth.json (used for submitting test jobs)

{
  "X_User_Email":"your email",
  "X_User_Token":"your user test token",
  "org_name": "your org name"
}

5. Run the protocol from the src dir

cd src
python protocols/bacteria_pelleting/protocol.py protocols/bacteria_pelleting/sample_input.json

## Publishing packages (requires unix)

You can add the packages to your transcriptic UI (allowing users to use them through the web UI).

1. Create a package by clicking "packages" and then "create new package" here https://secure.transcriptic.com/

2. Run the following commands from the build dir

python set_test_mode_false.py
bash make_zip.sh
transcriptic upload-release release.zip "YourPackageName"


 





