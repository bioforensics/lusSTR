# -------------------------------------------------------------------------------------------------
# Copyright (c) 2024, DHS.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------
#################################################################
#              Importing Necessary Packages                     #
#################################################################

import streamlit as st
from streamlit_option_menu import option_menu
import yaml
import subprocess
import os
import re

# ------ Packages For File/Folder Directory Selection --------- #

import tkinter as tk
from tkinter import filedialog

# Create a global Tkinter root window
root = tk.Tk()
root.withdraw()  # Hide the root window

#################################################################
#                        Functions                              #
#################################################################

# ------------ Function to Generate config.yaml File ---------- #
 
def generate_config_file(config_data, working_directory, workflow_type):
    if workflow_type == "STR":
        config_filename = 'config.yaml'
    elif workflow_type == "SNP":
        config_filename = 'snp_config.yaml'
    else:
        raise ValueError("Invalid workflow type. Please specify either 'STR' or 'SNP'.")

    config_path = os.path.join(working_directory, config_filename)
    with open(config_path, 'w') as file:
        yaml.dump(config_data, file)

# ------------ Function for folder selection ------------------ #

def folder_picker_dialog():
    folder_path = filedialog.askdirectory(master=root)
    return folder_path

# ------- Function for individual file selection -------------- #

def file_picker_dialog():
    file_path = filedialog.askopenfilename(master=root)
    return file_path

# ---- Function to validate prefix for output folder ---------- #

def validate_prefix(prefix):
    if re.match(r'^[A-Za-z0-9_-]+$', prefix):  # Allow alphanumeric characters, underscore, and hyphen
        return True
    else:
        return False

#################################################################
#              Front-End Logic For Navigation Bar               #
#################################################################

def main():
    
    # Page Layout (Theme and Fonts have been established in .streamlit/config.toml)
    st.set_page_config(layout='wide', initial_sidebar_state='collapsed')

    # Creating Navigation Bar
    
    selected = option_menu(
        menu_title=None,
        options=["Home", "STR", "SNP", "How to Use", "Contact"],
        icons=["house", "gear", "gear-fill", "book", "envelope"],
        menu_icon="cast",
        default_index=0,
        orientation="horizontal"
    )

    if selected == "Home":
        show_home_page()
        
    elif selected == "STR":
        show_STR_page()
        
    elif selected == "SNP":
        show_SNP_page()
        
    elif selected == "How to Use":
        show_how_to_use_page()
        
    elif selected == "Contact":
        show_contact_page()
        
#####################################################################
#                     lusSTR Home Page                              #
#####################################################################

def show_home_page():

    image_path = "logo.png"
    
    # CSS to hide full-screen button
    hide_img_fs = '''
    <style>
    button[title="View fullscreen"]{
    visibility: hidden;}
    </style>
    '''
    
    # Define column layout for centering image
    left_co, cent_co, last_co = st.columns([2.5, 8, 2.5])
    with cent_co:    
        st.image(image_path, use_column_width="auto")
    
    # Apply CSS to hide full-screen button
    st.markdown(hide_img_fs, unsafe_allow_html=True)

# -- Welcome Message Stuff

    st.markdown("""
        lusSTR is a tool written in Python to convert Next Generation Sequencing (NGS) data of forensic STR loci to different sequence representations (sequence bracketed form) and allele designations (CE allele, LUS/LUS+ alleles) for ease in downstream analyses.
        For more information on LusSTR, visit our [GitHub page](https://github.com/bioforensics/lusSTR/tree/master).
        """, unsafe_allow_html=True)

    st.info('Please Select One of the Tabs Above to Get Started on Processing Your Data!')

#####################################################################
#                        STR WORKFLOW                               #
#####################################################################

#####################################################################
# Specify STR Settings Which Will Be Used to Generate Config File   #
#####################################################################

def show_STR_page():

    st.title("STR Workflow")
    st.info('Please Select STR Settings Below for LusSTR! For Information Regarding the Settings, See the How to Use Tab.')
    
    # Input File Specification
    st.subheader("Input Files Selection")
    
    # Ask user if submitting a directory or individual file
    st.info("Please Indicate If You Are Providing An Individual Input File or a Directory Containing Multiple Input Files")
    input_option = st.radio("Select Input Option:", ("Individual File", "Directory with Multiple Files"))
    
    # Initialize session state if not already initialized
    if 'samp_input' not in st.session_state:
        st.session_state.samp_input = None
    
    # Logic for Path Picker based on user's input option

    if input_option == "Directory with Multiple Files":
        st.write('Please select a folder:')
        clicked = st.button('Folder Picker')
        if clicked:
            dirname = folder_picker_dialog()
            #st.text_input('You Selected The Following folder:', dirname)
            st.session_state.samp_input = dirname
            
    else:
        st.write('Please select a file:')
        clicked_file = st.button('File Picker')
        if clicked_file:
            filename = file_picker_dialog()
            #st.text_input('You Selected The Following file:', filename)
            st.session_state.samp_input = filename
    
    # Display The Selected Path
    if st.session_state.samp_input:
        st.text_input("Location Of Your Input File(s):", st.session_state.samp_input)
            
    # Store the Selected Path to Reference in Config
    samp_input = st.session_state.samp_input
                            
#####################################################################
#      STR: General Software Settings to Generate Config File       #
#####################################################################
   
    st.subheader("General Software")

    analysis_software = {'UAS': 'uas', 'STRait Razor v3': 'straitrazor', 'GeneMarker HTS': 'genemarker'}[st.selectbox("Analysis Software", options=["UAS", "STRait Razor v3", "GeneMarker HTS"], help="Indicate the analysis software used prior to lusSTR sex.")]

    sex = st.checkbox("Include sex-chromosome STRs", help = "Check the box if yes, otherwise leave unchecked.")
        
    output = st.text_input("Please Specify a prefix for generated output files or leave as default", "lusstr_output", help = "Be sure to see requirements in How to Use tab.")
                  
#####################################################################
#     STR: Convert Settings to Generate Config File                 #
#####################################################################
   
    st.subheader("Convert Settings")
     
    kit = {'ForenSeq Signature Prep': 'forenseq', 'PowerSeq 46GY': 'powerseq'}[st.selectbox("Library Preparation Kit", options=["ForenSeq Signature Prep", "PowerSeq 46GY"])]

    nocombine = st.checkbox("Do Not Combine Identical Sequences")
    
#####################################################################
#     STP: Filter Settings to Generate Config File                  #
#####################################################################
    
    st.subheader("Filter Settings")
    
    output_type = {'STRmix': 'strmix', 'EuroForMix': 'efm', 'MPSproto': 'mpsproto'}[st.selectbox("Output Type", options=["STRmix", "EuroForMix", "MPSproto"])]

    profile_type = {'Evidence': 'evidence', 'Reference': 'reference'}[st.selectbox("Profile Type", options=["Evidence", "Reference"])]

    data_type = {'Sequence': 'ngs', 'CE allele': 'ce', 'LUS+ allele': 'lusplus'}[st.selectbox("Data Type", options=["Sequence", "CE allele", "LUS+ allele"])]
 
    info = st.checkbox("Create Allele Information File")
    
    separate = st.checkbox("Create Separate Files for Samples", help = "If True, Will Create Individual Files for Samples; If False, Will Create One File with all Samples.")
    
    nofilters = st.checkbox("Skip all filtering steps", help = "Skip all Filtering Steps; Will Still Create EFM/MPSproto/STRmix Output Files")
    
    strand = {'UAS': 'uas', 'Forward': 'forward'}[st.selectbox("Strand Orientation", options=["UAS", "Forward"], help="Indicates the Strand Orientation in which to Report the Sequence in the Final Output Table for STRmix NGS only.")]
    
#####################################################################
#     STR: Specify Working Directory                                #
#####################################################################
    
    st.subheader("Set Working Directory")
    
    # Initialize session state if not already initialized
    if 'wd_dirname' not in st.session_state:
        st.session_state.wd_dirname = None
    
    clicked_wd = st.button('Please Specify A Working Directory Where You Would Like For All Output Results To Be Saved')
    if clicked_wd:
        wd = folder_picker_dialog()
        st.session_state.wd_dirname = wd
    
    # Display selected path
    if st.session_state.wd_dirname:
        st.text_input("Your Specified Working Directory:", st.session_state.wd_dirname)
    
    # Store Selected Path to Reference in Config
    wd_dirname = st.session_state.wd_dirname
    
#####################################################################
#     STR: Generate Config File Based on Settings                   #
#####################################################################
    
    # Submit Button Instance
    if st.button("Submit"):
    
        # Check if all required fields are filled
        if analysis_software and samp_input and output and wd_dirname:
    
            # Validate output prefix
            if not validate_prefix(output):
                st.warning("Please enter a valid output prefix. Only alphanumeric characters, underscore, and hyphen are allowed.")
                st.stop()  # Stop execution if prefix is invalid
    
            # Display loading spinner (Continuing Process Checks Above Were Passed)
            with st.spinner("Processing Your Data..."):
    
                # Construct config data
                
                config_data = {
                    "analysis_software": analysis_software,
                    "sex": sex,
                    "samp_input": samp_input,
                    "output": output,
                    "kit": kit,
                    "nocombine": nocombine,
                    "output_type": output_type,
                    "profile_type": profile_type,
                    "data_type": data_type,
                    "info": info,
                    "separate": separate,
                    "nofilters": nofilters,
                    "strand": strand
                }
    
                # Generate YAML config file
                generate_config_file(config_data, wd_dirname, "STR")
                
                # Subprocess lusSTR commands
                command = ["lusstr", "strs", "all"]
                
                # Specify WD to lusSTR
                if wd_dirname:
                    command.extend(["-w", wd_dirname + "/"])
    
                # Run lusSTR command in terminal
                try:
                    subprocess.run(command, check=True)
                    st.success("Config File Generated and lusSTR Executed Successfully! Output Files Have Been Saved to Your Designated Directory and Labeled with your Specified Prefix")
                except subprocess.CalledProcessError as e:
                    st.error(f"Error: {e}")
                    st.info("Please make sure to check the 'How to Use' tab for common error resolutions.")
    
        else:
            st.warning("Please make sure to fill out all required fields (Analysis Software, Input Directory or File, Prefix for Output, and Specification of Working Directory) before submitting.")
       
#####################################################################
#                        SNP WORKFLOW                               #
#####################################################################

#####################################################################
# Specify SNP Settings Which Will Be Used to Generate Config File   #
#####################################################################

def show_SNP_page():
    
    st.title("SNP Workflow")
    st.info('Please Select SNP Settings Below for lusSTR! For Information Regarding the Settings, See the How to Use Tab.')
    
    # Input File Specification
    st.subheader("Input Files Selection")
    
    # Ask user if submitting a directory or individual file
    st.info("Please Indicate If You Are Providing An Individual Input File or a Directory Containing Multiple Input Files")
    input_option = st.radio("Select Input Option:", ("Individual File", "Directory with Multiple Files"))
    
    # Initialize session state if not already initialized
    if 'samp_input' not in st.session_state:
        st.session_state.samp_input = None
    
    # Logic for Path Picker based on user's input option
    
    if input_option == "Directory with Multiple Files":
        st.write('Please select a folder:')
        clicked = st.button('Folder Picker')
        if clicked:
            dirname = folder_picker_dialog()
            #st.text_input('You Selected The Following folder:', dirname)
            st.session_state.samp_input = dirname
            
    else:
        st.write('Please select a file:')
        clicked_file = st.button('File Picker')
        if clicked_file:
            filename = file_picker_dialog()
            #st.text_input('You Selected The Following file:', filename)
            st.session_state.samp_input = filename
    
    # Display The Selected Path
    if st.session_state.samp_input:
        st.text_input("Location Of Your Input File(s):", st.session_state.samp_input)
    
    # Store Selected Path to Reference in Config
    samp_input = st.session_state.samp_input
            
#####################################################################
#      SNP: General Software Settings to Generate Config File       #
#####################################################################
    
    st.subheader("General Software")

    analysis_software = {'UAS': 'uas', 'STRait Razor v3': 'straitrazor'}[st.selectbox("Analysis Software", options=["UAS", "STRait Razor v3"], help="Indicate the analysis software used prior to lusSTR sex.")]

    output = st.text_input("Please Specify a prefix for generated output files or leave as default", "lusstr_output", help = "Be sure to see requirements in How to Use tab.")
     
    kit = {'Signature Prep': 'sigprep', 'Kintelligence': 'kintelligence'}[st.selectbox("Library Preparation Kit", options=["Signature Prep", "Kintelligence"])]

#####################################################################
#     SNP: Format Settings to Generate Config File                  #
#####################################################################
    
    st.subheader("Format Settings")
    
    # -- Select Type (Unique to SNP Workflow)
    types_mapping = {"Identify SNPs Only": "i", "Phenotype Only": "p", "Ancestry Only": "a", "All": "all"}
    selected_types = st.multiselect("Select Types:", options=types_mapping.keys(), help="Please Select a Choice or any Combination")
    types_string = "all" if "All" in selected_types else ", ".join(types_mapping.get(t, t) for t in selected_types)

    #if selected_types:
    #    st.text_input("You Selected:", types_string)
    
    # -- Filter
    nofilters = st.checkbox("Skip all filtering steps", help = "If no filtering is desired at the format step; if False, will remove any allele designated as Not Typed")

#####################################################################
#     SNP: Convert Settings to Generate Config File                 #
#####################################################################
    
    st.subheader("Convert Settings")
    
    separate = st.checkbox("Create Separate Files for Samples", help = "If want to separate samples into individual files for use in EFM")
    
    strand = {'UAS': 'uas', 'Forward': 'forward'}[st.selectbox("Strand Orientation", options=["UAS", "Forward"], help="Indicate which orientation to report the alleles for the SigPrep SNPs.")]
    
    # Analytical threshold value
    thresh = st.number_input("Analytical threshold value:", value=0.03, step=0.01, min_value = 0.0)

#####################################################################
#     SNP: Specify a Reference File if User Has One                 #
#####################################################################
    
    st.subheader("Specify a Reference File (Optional)")

    if 'reference' not in st.session_state:
        st.session_state.reference = None
   
    clicked_ref = st.button('Please Specify Your Reference File If You Have One', help = "List IDs of the samples to be run as references in EFM; default is no reference samples")
    if clicked_ref:
        ref = file_picker_dialog()
        st.session_state.reference = ref
    
    # Display Path to Selected Reference File
    if st.session_state.reference:
        st.text_input("Your Specified Reference File:", st.session_state.reference)
        
    # Store Selected Path to Reference in Config
    reference = st.session_state.reference
        
#####################################################################
#     SNP: Specify Working Directory                                #
#####################################################################
    
    st.subheader("Set Working Directory")    

    # Initialize session state if not already initialized
    if 'wd_dirname' not in st.session_state:
        st.session_state.wd_dirname = None
    
    clicked_wd = st.button('Please Specify A Working Directory Where You Would Like For All Output Results To Be Saved')
    if clicked_wd:
        wd = folder_picker_dialog()
        st.session_state.wd_dirname = wd
    
    # Display selected path
    if st.session_state.wd_dirname:
        st.text_input("Your Specified Working Directory:", st.session_state.wd_dirname)
    
    # Store Selected Path to Reference in Config
    wd_dirname = st.session_state.wd_dirname
    
#####################################################################
#     SNP: Generate Config File Based on Settings                   #
#####################################################################
    
    # Submit Button Instance
    if st.button("Submit"):
    
        # Check if all required fields are filled
        if analysis_software and samp_input and output and wd_dirname:
    
            # Validate output prefix
            if not validate_prefix(output):
                st.warning("Please enter a valid output prefix. Only alphanumeric characters, underscore, and hyphen are allowed.")
                st.stop()  # Stop execution if prefix is invalid
    
            # Display loading spinner (Continuing Process Checks Above Were Passed)
            with st.spinner("Processing Your Data..."):
    
                # Construct config data
                
                config_data = {
                    "analysis_software": analysis_software,
                    "samp_input": samp_input,
                    "output": output,
                    "kit": kit,
                    "types": types_string,
                    "thresh": thresh,
                    "separate": separate,
                    "nofilter": nofilters,
                    "strand": strand, 
                    "references": None  # Default value is None
                }
                
                # If a reference file was specified, add to config
                if reference:
                    config_data["references"] = reference
            
                # Generate YAML config file
                generate_config_file(config_data, wd_dirname, "SNP")
                
                # Subprocess lusSTR commands
                command = ["lusstr", "snps", "all"]
                
                # Specify WD to lusSTR
                if wd_dirname:
                    command.extend(["-w", wd_dirname + "/"])
    
                # Run lusSTR command in terminal
                try:
                    subprocess.run(command, check=True)
                    st.success("Config File Generated and lusSTR Executed Successfully! Output Files Have Been Saved to Your Designated Directory and Labeled with your Specified Prefix")
                except subprocess.CalledProcessError as e:
                    st.error(f"Error: {e}")
                    st.info("Please make sure to check the 'How to Use' tab for common error resolutions.")
    
        else:
            st.warning("Please make sure to fill out all required fields (Analysis Software, Input Directory or File, Prefix for Output, and Specification of Working Directory) before submitting.")
    
#####################################################################
#                        How To Use Page                            #
#####################################################################  
    
def show_how_to_use_page():
    
    st.title("Common Errors and Best Practices for Using lusSTR")
    
    st.header("1. File/Folder Path Formatting")
   
    st.write("Please ensure that the displayed path accurately reflects your selection. When using the file or folder picker, navigate to the desired location and click 'OK' to confirm your selection.")

    st.header("2. Specifying Output Prefix")
    
    st.write("The purpose of specifying the output prefix is for lusSTR to create result files and folders with that prefix in your working directory. Please ensure that you are following proper file naming formatting and rules when specifying this prefix. Avoid using characters such as '/', '', '.', and others. Note: To avoid potential errors, you can simply use the default placeholder for output.")
    
    st.code("Incorrect: 'working_directory/subfolder/subfolder'\nCorrect: working_directory/output # or just output, since you will likely already be in the working directory when specifying it before submitting.")
    
    st.write("Note that some result files may be saved directly in the working directory with the specified prefix, while others will be populated in a folder labeled with the prefix in your working directory.")
    st.write("Be aware of this behavior when checking for output files.")

    st.header("3. Specifying Working Directory")
    st.write("Please Ensure That You Properly Specify a Working Directory. This is where all lusSTR output files will be saved. To avoid potential errors, specifying a working directory is required.")
    
    st.title("About lusSTR")
    
    st.markdown("""
                        
    **_lusSTR Accommodates Four Different Input Formats:_**

    (1) UAS Sample Details Report, UAS Sample Report, and UAS Phenotype Report (for SNP processing) in .xlsx format (a single file or directory containing multiple files)

    (2) STRait Razor v3 output with one sample per file (a single file or directory containing multiple files)

    (3) GeneMarker v2.6 output (a single file or directory containing multiple files)

    (4) Sample(s) sequences in CSV format; first four columns must be Locus, NumReads, Sequence, SampleID; Optional last two columns can be Project and Analysis IDs.

   
    """, unsafe_allow_html = True)
    
#####################################################################
#                        Contact Page                            #
#####################################################################  

def show_contact_page():
    st.title("Contact Us")
    st.write("For any questions or issues, please contact rebecca.mitchell@st.dhs.gov, daniel.standage@st.dhs.gov, or s.h.syed@email.msmary.edu")

#####################################################################
#     lusSTR RUN                                                    #
#####################################################################

if __name__ == "__main__":
    main()

def subparser(subparsers):
    parser = subparsers.add_parser("gui", help="Launch the Streamlit GUI")
    parser.set_defaults(func=main)
