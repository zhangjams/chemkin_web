import tkinter
from tkinter import messagebox


def helloCallBack():
   messagebox.showinfo( "Hello Python", "Hello World")

def start_gui():
	top = tkinter.Tk()
	B = tkinter.Button(top, text ="Hello", command = helloCallBack)
	B.pack()

	top.mainloop()