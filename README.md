# delve_tx_public

This is a set of transcriptic package protocols.

## Protocol Overview



## Running protocols from the command line


1. update protocols/<the_protocol>/sample_input.json with the correct container id's you want to use
2. update ~/.transcriptic to have your login credentials

{
  "api_root": "https://secure.transcriptic.com",
  "token": "<your token>",
  "email": "<your email>",
  "organization_id": "<your org>"
}

3. Add the following to src/auth.json

{
  "X_User_Email":"<your email>",
  "X_User_Token":"<your user token>",
  "org_name": "<your org name>"
}

4. Add the following to src/test_mode_auth.json (used for submitting test jobs)

{
  "X_User_Email":"<your email>",
  "X_User_Token":"<your user test token>",
  "org_name": "<your org name>"
}

5. Run the protocol

cd src
python protocols/



## Publishing packages


 





