__author__ = 'Paulina Kania'

from Tkinter import *
import tkFileDialog
import tkMessageBox
from PIL import Image
from PIL import ImageTk
import sys
import os

import HMM
import input
from ChouFasman import ChouFasman
import comparator
import DSSP_nasz


import BiGRU

app = Tk()
app.title("Secondary Structure Prediction")
app.geometry('800x650+200+200')
app.configure(background="beige")


# options definition for opening or saving a file
app.file_opt_in = options = {}
options['defaultextension'] = '.fasta'
options['filetypes'] = [('FASTA files', '.fasta'), ('Output', '.out'), ('All files', '.*')]
options['initialdir'] = 'C:\\'
options['initialfile'] = 'myfile.fasta'
options['parent'] = app
options['title'] = 'Open file'

app.file_opt_out = options = {}
options['defaultextension'] = '.fasta'
options['filetypes'] = [('Output', '.output')]
options['initialdir'] = 'C:\\'
options['initialfile'] = 'out.output'
options['parent'] = app
options['title'] = 'Save file'

def restart_program():
    """Restarts the current program.
    Note: this function does not return. Any cleanup action (like
    saving data) must be done before calling this function."""
    python = sys.executable
    os.execl(python, python, * sys.argv)

def select_all(event):
    """
    This function make the selection of the inserted text into the text field
    :param event
    """
    text.tag_add(SEL, "1.0", END)
    text.mark_set(INSERT, "1.0")
    text.see(INSERT)
    return 'break'

def clear_text():
    """
    This function clears the entry in the text field
    """
    text.delete(1.0, END)

def browseFile():
    """
    This function lets the user choose the file to load;
    after that the loaded file is put into the variable 'text',
    which represents the white window in the GUI
    :return: the input in FASTA format
    """
    app.withdraw() #hides the browse window
    file_path = tkFileDialog.askopenfile(mode='r', **app.file_opt_in)
    file = file_path.read()
    app.deiconify() #wraca do glownego okna, w stanie w jakim zostalo przed wczytaniem pliku
    text.insert(END, file) #wpisuje tekst z pliku do okna

def save_file():
     """Returns an opened file in write mode."""
     app.withdraw() #hides the browse window
     file_path = tkFileDialog.asksaveasfile(mode='w', **app.file_opt_out)
     app.deiconify() #wraca do glownego okna, w stanie w jakim zostalo przed wczytaniem pliku
     return file_path

def help_index():
    """
    This function shows the window with help for the user
    """
    tkMessageBox.showinfo("Help", "This is a message of help")

def show_about():
    """
    This function shows some information about this program, for example about authors
    """
    tkMessageBox.showinfo("About", "Secondary Structure Prediction program")

menubar = Menu(app)
menubar.configure(background="lightgreen")
filemenu = Menu(menubar, tearoff=0)
filemenu.add_command(label="New", command=restart_program)
filemenu.add_command(label="Open", command=browseFile)
filemenu.add_command(label="Save as...", command=save_file)

filemenu.add_separator()

filemenu.add_command(label="Exit", command=app.quit)
menubar.add_cascade(label="File", menu=filemenu)

helpmenu = Menu(menubar, tearoff=0)
helpmenu.add_command(label="Help", command=help_index)
helpmenu.add_command(label="About...", command=show_about)
menubar.add_cascade(label="Help", menu=helpmenu)

app.config(menu=menubar)

labelText = StringVar()
labelText.set("Insert FASTA code here")
label1 = Label(app, textvariable=labelText, height=2)
label1.grid(row=0, column=0, columnspan=2, sticky='W', padx=10, pady=5)
label1.configure(background="beige")
button1 = Button(app, text='Browse', width=10, command=browseFile)
button1.configure(background="lightgreen")
button1.grid(row=0, column=1, sticky='W', pady=5)

# app.entryVariable = StringVar()
# app.E1 = Entry(app,textvariable=app.entryVariable, bd=5,width=17)
# app.E1.configure(background="lightgreen")
# app.E1.grid(row=0, column=2, sticky='W')
# #root.E1.bind("<Return>",root.OnPressEnter)
# app.entryVariable.set('Enter PDB ID here')
# app.E1.focus_set()
# app.E1.selection_range(0, END)

text = Text(app, wrap=WORD, width=70, height=20)
text.grid(row=1, column=0, columnspan=3, padx=10, sticky='W')
text.bind("<Control-Key-a>", select_all)
text.bind("<Control-Key-A>", select_all) # just in case caps lock is on

scrollbar = Scrollbar(app)
yscroll = Scrollbar(command=text.yview, orient=VERTICAL)
yscroll.configure(background="lightgreen")
yscroll.grid(row=1, column=2, sticky='NSE')
text.configure(yscrollcommand=yscroll.set)

labelText2 = StringVar()
labelText2.set("Methods")
label2 = Label(app, textvariable=labelText2, height=2)
label2.configure(background="beige")
label2.grid(row=2, column=0, sticky='WN', padx=10)
CheckVar1 = IntVar()
CheckVar2 = IntVar()
CheckVar3 = IntVar()
CheckVar4 = IntVar()
C1 = Checkbutton(app, text = "HMM", variable = CheckVar1,onvalue = 1, offvalue = 0, height=3, width = 15)
C2 = Checkbutton(app, text = "Statistical", variable = CheckVar2,onvalue = 1, offvalue = 0, height=3, width = 15)
C3 = Checkbutton(app, text = "Consensus", variable = CheckVar3,onvalue = 1, offvalue = 0, height=3, width = 15)
C4 = Checkbutton(app, text = "BiGRU", variable = CheckVar4,onvalue = 1, offvalue = 0, height=3, width = 15)
C1.grid(row=3, column=0, sticky='N')
C2.grid(row=3, column=1, sticky='N')
C3.grid(row=3, column=2, sticky='N')
C4.grid(row=3, column=3, sticky='N')


labelText3 = StringVar()
labelText3.set("Output options")
label3 = Label(app, textvariable=labelText3)
label3.configure(background="beige")
label3.grid(row=4, column=0, sticky='WN', padx=10)
CheckVar4 = IntVar()
CheckVar5 = IntVar()
CheckVar6 = IntVar()
CheckVar7 = IntVar()
CheckVar8 = IntVar()
CheckVar9 = IntVar()
CheckVar10 = IntVar()
C4 = Checkbutton(app, text = "Aminoacids sequence", variable = CheckVar4, onvalue = 1, offvalue = 0, height=3, width = 20)
C5 = Checkbutton(app, text = "Chou-Fasman", variable = CheckVar5, onvalue = 1, offvalue = 0, height=3, width = 20)
C6 = Checkbutton(app, text = "HMM Structure", variable = CheckVar6, onvalue = 1, offvalue = 0, height=3, width = 20)
C7 = Checkbutton(app, text = "DSSP Structure", variable = CheckVar7, onvalue = 1, offvalue = 0, height=3, width = 20)
C8 = Checkbutton(app, text = "Consensus Structure", variable = CheckVar8, onvalue = 1, offvalue = 0, height=3, width = 20)
C9 = Checkbutton(app, text = "BiGRU Structure", variable = CheckVar9, onvalue = 1, offvalue = 0, height=3, width = 20)
C10 = Checkbutton(app, text = "Scoring", variable = CheckVar10, onvalue = 1, offvalue = 0, height=3, width = 20)

C4.grid(column=0, row=5, sticky='W')
C5.grid(column=0, row=6, sticky='W')
C6.grid(column=0, row=7, sticky='W')
C7.grid(column=1, row=5, sticky='W')
C8.grid(column=1, row=6, sticky='W')
C9.grid(column=0, row=8, sticky='W')
C10.grid(column=1, row=7, sticky='W')

C1.configure(background="beige", activebackground="lightgreen")
C2.configure(background="beige", activebackground="lightgreen")
C3.configure(background="beige", activebackground="lightgreen")
C4.configure(background="beige", activebackground="lightgreen")
C5.configure(background="beige", activebackground="lightgreen")
C6.configure(background="beige", activebackground="lightgreen")
C7.configure(background="beige", activebackground="lightgreen")
C8.configure(background="beige", activebackground="lightgreen")
C9.configure(background="beige", activebackground="lightgreen")
C10.configure(background="beige", activebackground="lightgreen")
def appRUN():
    """
    this function starts when the user presses the "Run" button;
    it creates the output file according to the options and the methods chosen by the user
    """
    fastaCode = text.get(1.0, END) #pobiera tresc znajdujaca sie w oknie (kod FASTA)

    pdbID, seqInput = input.check_input(fastaCode)
    CF = ChouFasman()
    GRU = GRUNetwork()
    #ss_viterbi = ""
    ss_forward = ""
    ss_backward = ""
    ss_chou = ""
    consensus = ""
    #pdbID = app.E1.get()
    ss_dssp = DSSP_nasz.dssp_ss(pdbID)
    ss_dssp_out = ""
    seqInput_out = ""
    diff = len(ss_dssp) - len(seqInput)
    ss_BiGRU = ""
    # CHECKING LENGTH OF INPUT
    if diff > 0:
        for i in range(len(seqInput)):
            ss_dssp_out += ss_dssp[i]
            seqInput_out += seqInput[i]
    elif diff < 0:
        result = tkMessageBox.askokcancel("Warning", "Your input sequence is longer than the"
                                         "DSSP structure in our database."
                                         "Do you agree to cut your sequence? If you click cancel "
                                         "the predicted sequences may have different length")
        if result is True:
            for i in range(len(ss_dssp)-1):
                ss_dssp_out += ss_dssp[i]
                seqInput_out += seqInput[i]
        if result is False:
            ss_dssp_out = ss_dssp
            seqInput_out = seqInput
    else:
        ss_dssp_out = ss_dssp
        seqInput_out = seqInput
    out_file = save_file()
    if CheckVar4.get() == 1:
        out_file.write("Aminoacids sequence \n" + seqInput_out)
    if CheckVar1.get() == 1 and CheckVar6.get() == 1:
        #exec HMM and write the sequences into the output file
        ss_forward = HMM.mapping_num_to_letters(HMM.hmm_structure('forward', seqInput_out))
        ss_backward = HMM.mapping_num_to_letters(HMM.hmm_structure('backward', seqInput_out))
        out_file.write("\n\n--- HMM Structure --- Forward algorithm ---\n" + ss_forward)
        out_file.write("\n\n--- HMM Structure --- Backward algorithm ---\n" + ss_backward)
    if CheckVar2.get() == 1 and CheckVar5.get() == 1:
        # exec Chou-Fasman and write the sequences into the output file
        ss_chou = ChouFasman.engine(CF, seq_fasta=seqInput_out)
        out_file.write("\n\n--- Chou-Fasman ---\n" + ss_chou)
    if CheckVar2.get() == 1:
        #if Chou Fasman method is selected
        ss_chou = ChouFasman.engine(CF, seq_fasta=seqInput_out)
    if CheckVar1.get() == 1:
        #if HMM method is selected
        ss_forward = HMM.mapping_num_to_letters(HMM.hmm_structure('forward', seqInput_out))
        ss_backward = HMM.mapping_num_to_letters(HMM.hmm_structure('backward', seqInput_out))
    if CheckVar3.get() == 1 and CheckVar8.get() == 1:
        #if Consensus and Consensus Structure fields are selected
        ss_forward = HMM.mapping_num_to_letters(HMM.hmm_structure('forward', seqInput_out))
        ss_backward = HMM.mapping_num_to_letters(HMM.hmm_structure('backward', seqInput_out))
        ss_chou = ChouFasman.engine(CF, seq_fasta=seqInput_out)
        list_of_seq=[ss_forward,ss_backward,ss_chou]
        if len(list_of_seq) != 3:
            tkMessageBox.showwarning("Warning", "Please check if you have selected all the methods")
        else:
            consensus = comparator.HMM_consensus(list_of_seq)
            out_file.write("\n\n--- Consensus structure ---\n" + consensus)
    if CheckVar7.get() == 1:
        out_file.write("\n\n--- DSSP Structure ---\n" + ss_dssp_out)
    if CheckVar9.get() == 1:
        ss_BiGRU = GRU.secondary_structure(seqInput_out)
        out_file.write("\n\n--- BiGRU structure ---\n" + ss_BiGRU)
    if CheckVar10.get() == 1:
        ss_forward = HMM.mapping_num_to_letters(HMM.hmm_structure('forward', seqInput_out))
        ss_backward = HMM.mapping_num_to_letters(HMM.hmm_structure('backward', seqInput_out))
        ss_chou = ChouFasman.engine(CF, seq_fasta=seqInput_out)
        consensus = comparator.HMM_consensus([ss_forward,ss_backward,ss_chou])
        score1 = comparator.identities([ss_forward, ss_dssp_out])
        score2 = comparator.identities([ss_backward, ss_dssp_out])
        score3 = comparator.identities([ss_chou, ss_dssp_out])
        score4 = comparator.identities([consensus, ss_dssp_out])
        score1 = str(score1)
        score2 = str(score2)
        score3 = str(score3)
        score4 = str(score4)
        out_file.write("\n\n\n******* IDENTITIES PERCENTAGE *******\n\n" + "HMM Forward to DSSP:\t" + score1 + \
                        "\nHMM Backward to DSSP:\t" + score2 + "\nChou-Fasman to DSSP:\t" + score3 + \
                        "\nConsensus to DSSP:\t" + score4)
        percentage = comparator.compare_structures([ss_chou, ss_backward, ss_forward, ss_dssp_out, consensus])
        out_file.write("\n\n\n******* STRUCTURE TYPE PERCENTAGE *******\n\n" \
                       + '\t         H\t    E\t    T\t  C\n' \
                       + "Chou-Fasman: \t" + str(percentage[0])
                       + "\nHMM Backward:\t" + str(percentage[1])
                       + "\nHMM Forward: \t" + str(percentage[2])
                       + "\nDSSP:        \t" + str(percentage[3])
                       + "\nConsensus:   \t" + str(percentage[4]))

button2 = Button(app, text='Run', width=10, command=appRUN)
button2.configure(background="lightgreen")
button2.grid(row=8, column=2, sticky='S')

button5 = Button(app, text='Reset', width=10, command=clear_text)
button5.configure(background="lightgreen")
button5.grid(row=7, column=2, sticky='S')

def open_output():
    app.withdraw() #hides the browse window
    file_path = tkFileDialog.askopenfile(mode='r', **app.file_opt_out)
    file = file_path.read()
    app.deiconify() #wraca do glownego okna, w stanie w jakim zostalo przed wczytaniem pliku
    return text.insert(END, file)

button4 = Button(app, text='Open output', width=10, command=open_output)
button4.configure(background="lightgreen")
button4.grid(row=7, column=3, sticky='S')

button4 = Button(app, text='Exit', width=10, command=app.quit)
button4.configure(background="lightgreen")
button4.grid(row=8, column=3, sticky='S')

img = ImageTk.PhotoImage(Image.open("logoSSP6.jpg"))
imglabel = Label(app, image=img).grid(row=0, column=3, rowspan=3)


labelText3 = StringVar()
labelText3.set("Choose or paste a file in FASTA format; \n"
               "then choose the methods "
               "you want \nto use and the options for the output. \nThen run the application.\n"
               "Use only unambiguous protein \nalphabet! "
               "You will receive \nan empty output. Scoring \n"
               "does not work without PDB ID.")
label3 = Label(app, textvariable=labelText3, width=30, wraplength=350, height=10, anchor=W, justify=LEFT)
label3.configure(background="beige")
label3.grid(column=3, row=3, rowspan=5, sticky='N')

if __name__ == '__main__':
    app.mainloop()
