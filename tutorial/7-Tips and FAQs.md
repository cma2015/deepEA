<div align='center' ><font size='70'>Tips and FAQs</font></div>

## Tips

The folloing screenshot shows us:

- **How to view the data in deepEA**
- **How to re-run a job**
- **How to save data to your local disk**

	![0-0](../assets/img/0-0.png)



## Frequently Asked Questions (FAQs)

### How to upload in-house data into deepEA local server

<a href="https://youtu.be/vDd9yQHiYYQ" target="_blank">
    <img border="0" src="../assets/img/how_to_upload_data.png" />
</a>

### How to become an admin user

- First, register with email `admin@example.org`, set the password arbitrary.

- Then, login with `admin@example.org`


### How to stop deepEA local server
- Press `Ctrl + C` (for windows and unix users) or `Cmd + C` (for Mac OS users)

### How to re-launch deepEA local server when I exited the docker container
- First, using the following command to check the container ID
  ```bash
  docker ps -a
  ```
- Then run the following command
  ```bash
  docker container start  container ID
  docker exec -it container ID bash
  bash /home/galaxy/run.sh
  ```
### How to mount local disk into deepEA docker container

```bash
docker run -it -v /your home directory:/home/galaxy/database/files/000 -p 8080:8080 malab/deepea bash
```



