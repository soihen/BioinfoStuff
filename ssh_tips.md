## simplify ssh
generate a so-called `config` file in the folder ~/.ssh 

The config should look like this:
```
Host kai
  Hostname 123.456.78.910
  User kai
  ForwardX11 yes
```
Then change the mode:
```
[kai@admin]$ chmod 700 config
```

---

## skip password
1. Create public and private keys using ssh-key-gen on local-host
```
[kai@admin]$ ssh-keygen -t rsa
```
then follow the instruction printed on screen

You could press the enter key directly without having a passphrase

2. Copy the public key to remote-host using ssh-copy-id
```
[kai@admin]$ ssh-copy-id -i ~/.ssh/id_rsa.pub remote-host
```
enter the password on remote-host
