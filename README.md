# cfht2caom2
Starting point to build an application to generate CAOM2 Observations from FITS files.

# How to Run CFHT

In an empty directory (the 'working directory'), on a machine with Docker installed:

1. In the master branch of this repository, find the scripts directory, and copy the file cfht_run_state.sh to the working directory. e.g.:

  ```
  wget https://raw.github.com/opencadc-metadata-curation/cfht2caom2/master/scripts/cfht_run_state.sh
  ```

2. Ensure the script is executable:

```
chmod +x cfht_run_state.sh
```

3. To run the application:

```
./cfht_run_state.sh
```

