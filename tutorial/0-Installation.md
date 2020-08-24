## deepEA installation

- **Step 1**: Docker installation

  **i) Docker installation and start (<a href="https://docs.docker.com/install" target="_blank">Official installation tutorial</a>)**

  For **Windows (Only available for Windows 10 Prefessional and Enterprise version):**

	* Download <a href="https://download.docker.com/win/stable/Docker%20for%20Windows%20Installer.exe" target="_blank">Docker</a> for windows;
	* Double click the EXE file to open it;
  * Follow the wizard instruction and complete installation;
  * Search docker, select ___Docker for Windows___ in the search results and click it.
  
  For **Mac OS X (Test on macOS Sierra version 10.12.6 and macOS High Sierra version 10.13.3):**
  
	- Download <a href="https://download.docker.com/mac/stable/Docker.dmg" target="_blank">Docker</a> for Mac OS;
  * Double click the DMG file to open it;
  * Drag the docker into Applications and complete installation;
  * Start docker from Launchpad by click it.

	For **Ubuntu (Test on Ubuntu 18.04 LTS):**

  * Go to <a href="https://download.docker.com/linux/ubuntu/dists/" target="_blank">Docker</a>, choose your Ubuntu version, browse to **pool/stable** and choose **amd64, armhf, ppc64el or s390x**. Download the **DEB** file for the Docker version you want to install;
  * Install Docker, supposing that the DEB file is download into following path:___"/home/docker-ce<version-XXX>~ubuntu_amd64.deb"___ </br>

    ```bash
      $ sudo dpkg -i /home/docker-ce<version-XXX>~ubuntu_amd64.deb      
      $ sudo apt-get install -f
    ```


  **ii) Verify if Docker is installed correctly**

  Once Docker installation is completed, we can run `hello-world` image to verify if Docker is installed correctly. Open terminal in Mac OS X and Linux operating system and open CMD for Windows operating system, then type the following command:

  ```bash
 $ docker run hello-world
  ```

   **<font color =red>Note</font>:** root permission is required for Linux operating system.

- Once Docker is installed successfully, you will see the following message:
  

![docker](../assets/img/Verify_docker.png)
	
- **Step 2**: deepEA installation from Docker Hub
```bash
# pull latest deepEA Docker image from docker hub
$ docker pull malab/deepea
```
- **Step 3**: Launch deepEA local server
```bash
$ docker run -it -p 8080:8080 malab/deepea bash
$ bash /home/galaxy/run.sh
```

Then, deepea local server can be accessed via http://localhost:8080

![index](../assets/img/index.png)

