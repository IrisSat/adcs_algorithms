# adcs_algorithms
Windows Users:
1. Required programs to build the project:
    - 1.1 Install Visual Studio Code: https://code.visualstudio.com/download
    - 1.2 Install MSYS2: https://www.msys2.org/
    - 1.3 Install make: https://sourceforge.net/projects/gnuwin32/
        - http://gnuwin32.sourceforge.net/packages/make.htm
    - 1.4 Add the following to your Path environment variable: C:\Program Files (x86)\GniWin32\bin
        - https://docs.microsoft.com/en-us/previous-versions/office/developer/sharepoint-2010/ee537574(v=office.14)
2. Install the following Visual Studio Code Extensions: https://code.visualstudio.com/docs/editor/extension-marketplace
    - 2.1 C/C++
    - 2.1 Makefile Tools
3. Set up your project for execution/debugging:
    - 3.1 Open Makefile Tools extenstion (press "Makefile" icon, left column, below explorer, search, run and debug...)
        **Note: you may have to press F5 on first try (with ADCS.c file open/active) or ctrl+shift+D or "Run and Debug" icon for the Makefile icon to appear (the build will fail)
    - 3.2 Click "..." icon --> "Makefile: Configure"
    - 3.3 Set "Build Target" to "all"
    - 3.4 Click "Makefile: Build the current target" icon
    - 3.5 Set "Launch Target" to {project path}\adcs_algorithms?ADCS()
    - 3.6 Press Debug or Play icon to execute code
