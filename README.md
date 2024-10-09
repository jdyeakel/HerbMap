# Diffision mapping of herbivore and primate morphology Julia Package

To develop locally, clone repo and set up julia package in development mode:

Yes, you are correct again! After initiating `dev` mode, you need to navigate to the package directory using `cd(path_to_project)` before running `instantiate`. Here's how it works step-by-step:

### **Correct Workflow:**

1. **Clone the Repository**:

   For SSH:
   ```bash
   git clone git@github.com:jdyeakel/HerbMap.git
   ```

   For HTTPS:
   ```bash
   git clone https://github.com/jdyeakel/HerbMap.git
   ```

2. **Open Julia**:

   In the terminal, start Julia:

   ```bash
   julia
   ```

3. **Activate `dev` Mode**:

   Press `]` to enter package mode and activate the `dev` mode by specifying the full path to the `HerbMap` project directory:

   ```julia
   ]
   dev /path/to/HerbMap
   ```

4. **Navigate to the Project Directory**:

   After activating `dev` mode, return to the Julia prompt by pressing **backspace** (to exit package mode) and then navigate to the project directory:

   ```julia
   cd("/path/to/HerbMap")
   ```

5. **Instantiate the Environment**:

   Once you are in the project directory, instantiate the environment to install the dependencies:

   ```julia
   pkg> instantiate
   ```

### **Full Correct Workflow for Your Collaborator**

1. **Clone the repository**:
   ```bash
   git clone git@github.com:jdyeakel/HerbMap.git
   ```

2. **Open Julia** and activate `dev` mode:
   ```julia
   julia
   ]
   dev /path/to/HerbMap  # Activate dev mode
   ```

3. **Navigate to the project directory**:
   ```julia
   cd("/path/to/HerbMap")  # Change directory to the project folder
   ```

4. **Instantiate the environment**:
   ```julia
   pkg> instantiate  # Install dependencies
   ```

This ensures everything is properly set up for development. Thanks for pointing that out! Let me know if you need anything else.
