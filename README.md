# Getting Started

`fdaPDE` is a C++ library that interfaces with **R**, one of the most popular languages for data analysis. It runs on all major platforms: **Linux**, **macOS**, and **Windows**.

The development version of the library has **strict system requirements**, which may be difficult to satisfy manually. In particular, the package must be built using **GCC version 15**, which is not available in some commonly used Linux distributions, such as Ubuntu, nor is it easily installed on macOS. To simplify installation, we provide a Docker image containing a pre-configured environment running the R package [`fdaPDE-R`](https://github.com/fdaPDE/fdaPDE-R).

---

## Recommended Setup

* **macOS** and **Linux** users: We strongly recommend using the Docker image.
* **Windows** users: You can either:

  * If your machine supports **R version ≥ 4.5.0**, install it along with [RTools45](https://cran.r-project.org/bin/windows/Rtools/). This is the preferred approach, if compatible.
  * Alternatively, use the Docker image (requires WSL – see instructions below).

---

## Image Usage

If you have not already installed Docker on your machine, please refer to the **Installation** section below.
To pull the Docker image, run the following command from a terminal:

```bash
docker pull aldoclemente/fdapde-docker:rstudio-nographics
```

To run a container, execute:

```bash
docker run --rm -d -p 8787:8787 -v /path/to/data:/home/user/data --name rstudio -e PASSWORD=password aldoclemente/fdapde-docker:rstudio-nographics
```

* Replace `/path/to/data` with the **full path** to the folder containing your data or course material.

This will launch an **RStudio Server** instance inside Docker. You can then access RStudio in your browser at:

**[http://localhost:8787](http://localhost:8787)**

Log in with:

* **Username**: `user`
* **Password**: `password`

Inside RStudio, set your working directory by running:

```
setwd("data/")
```

You have read/write permissions in this directory, so you can modify scripts, save results, etc. All changes are reflected on your local machine. You can safely close and reconnect to [http://localhost:8787](http://localhost:8787) at any time.

To stop the container when you're done:

```
docker stop rstudio
```

---

## Docker Installation

### Windows Users

#### Step 1: Enable Virtualization

1. Restart your computer and enter the **BIOS setup** (usually via F2, F10, Del, or Esc).
2. Enable **virtualization** (e.g., Intel VT-x or AMD-V).
3. Save and exit the BIOS.

#### Step 2: Enable WSL and Virtual Machine Platform

1. Open the **"Turn Windows features on or off"** menu (search via Start or Cortana).
2. Enable:

   * **Windows Subsystem for Linux**
   * **Virtual Machine Platform**

#### Step 3: Install Ubuntu via WSL

1. Open **Command Prompt** or **PowerShell** as Administrator.
2. Run:

```
wsl --install --distribution Ubuntu
```

3. Restart your machine and follow the instructions to set up your Ubuntu user account.
4. Verify your WSL version:

```
wsl --version
```

#### Step 4: Install Docker Desktop

1. Download and install Docker Desktop from [docker.com](https://www.docker.com/).
2. Log in or create a Docker Hub account.
3. Go to **Settings > Resources > WSL Integration** and enable integration for Ubuntu.
4. Restart your machine.

#### Step 5: Test Docker

1. Open your WSL terminal by running from PowerShell:

```
wsl
```

2. Run:

```
docker --version
```

---

### macOS Users

1. Download and install Docker Desktop from [this link](https://docs.docker.com/desktop/setup/install/mac-install/).
2. Open Docker Desktop and grant the necessary permissions.
3. Log in or create a Docker Hub account.
4. Verify the Docker installation:

```bash
docker --version
```

---

### Ubuntu Users

1. Install Docker:

```
sudo apt-get update
sudo apt-get install -y docker.io
```

2. Check that Docker is running:

```
systemctl status docker
```

If Docker is inactive, start the service:

```
sudo systemctl start docker
```

3. To avoid typing `sudo` whenever you run the `docker` command, add your username to the Docker group:

```
sudo usermod -aG docker $USER
```

4. Restart your machine.

5. Log in to Docker Hub:

```
docker login
```

6. You can now pull the Docker image and run the container as described above.

---

## Final Notes

* Ensure Docker is running before pulling or running any image.
* If you encounter issues, consult the [official Docker documentation](https://docs.docker.com/) for your platform.


