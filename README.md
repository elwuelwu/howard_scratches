# Mathematical scratches for decoding binary chirps

The document matrices.tex contains matrix language research of binary chirps. There are also some Julia and Mathematica code for practical simulation purposes.

## Workflow

Clone the reprository by

```bash
git clone git@github.com:elwuelwu/howard_scratches.git
```

Create your own working branch, push the changes to remote. When the added feature is ready, rebase the branch from master and create a pull request.

## Running the code in Aalto computing services

### Opening the virtual desktop

Search for VMware Horizon Client, connect to vdi.aalto.fi and launch Ubuntu session. 

Do the following steps when you have opened the Ubuntu session **only for the first time**, these steps need only be ran once:

1. SSH key generation: In the Ubuntu session, open terminal with `CTRL+ALT+T` and type `ssh-keygen` Press enter to all queries so that the ssh key is stored in the default location. Just press enter when asked for passphrase so no passwords are set.
2. SSH key authorizing: You are connecting between Aalto machines and all the machines are using the same keys. Thus you need to authenticate the access with your own public key.
This is done by first typing in the terminal `cd ~/.ssh` then `cat id_rsa.pub`. Copy the printed SSH key and type in `nano authorized_keys` to open a file called "authorized_keys" in the nano text editor. Paste the copied key and save by `CTRL+X` and type `Y` when asked if you want to save the file.

Now we have a seamless SSH access between the Aalto computers.

### Connecting to the Brute & Force computing servers

Aalto has two servers, namely Brute and Force, where one can run long simulations for free. The same instructions work for both servers. The address for Brute is `brute.aalto.fi` and for Force it is `force.aalto.fi`.
I am going to use Brute as an example since it usually has less load.

Type in `ssh -X brute.aalto.fi` in order to connect to Brute server with the X11 forwarding that allows you to remotely use GUI applications. There might be a query whether you trust the computer you are connecting to. Continue connecting and it should be able to connect seamlessly with the generated SSH keys.

Now we can start our favorite GUI applications. Type in for example
`gnome-system-monitor &` (be patient, this command might take a while to complete)
The `&`-sign determines that the process is ran on the background of the terminal, so that the terminal is not occupied.
From the system monitor you can easily see the load of the server at the moment. If the server is overly loaded, try the above steps with swapping to force server instead of brute.
Now we can start running the simulations with
`mathematica &`
or
`matlab &`
ect.

You can also detach from the session easily by simply closing the VMware window and even close your local computer. The session stays alive in the cloud for 24 hours, so in the case of longer simulations, you should at least once a day connect to the session.
You can re-attach to the session simply by opening VMware client and launching the same Ubuntu session again.
