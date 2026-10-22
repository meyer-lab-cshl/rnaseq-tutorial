# Running the tutorial in Docker

The Docker image contains snakemake and every conda environment in `envs/`,
so you do not need to install miniforge or snakemake. Install
[Docker Desktop](https://www.docker.com/products/docker-desktop/) and give it
at least 4 GB of memory in Settings > Resources.

## Command line

Build the image once, from inside this folder:

```bash
docker build -t rnaseq-tutorial .
```

Start a shell in the container, with this folder mounted at `/work`:

```bash
docker run --rm -it -v "$PWD":/work rnaseq-tutorial
```

Test the setup inside the container:

```bash
snakemake --cores 1 -s setup.smk
```

## VS Code

Install the
[Dev Containers](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers)
extension, open this folder in VS Code, and choose "Reopen in Container". The
first start builds the image. Terminals in VS Code then run inside the
container, where you can run the same test:

```bash
snakemake --cores 1 -s setup.smk
```
