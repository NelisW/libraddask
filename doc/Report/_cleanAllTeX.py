#clean a precisely defined set of files from a precisely defined directory set.
#the specified files are deleted, after the user has been prompted


import os.path, fnmatch
import sys
from optparse import OptionParser
import pyradi.ryfiles as ryfiles


#function to delete files in specified directory
#  first parameter defines if the search must be recursively 0=not, 1=recursive
#  second parameter specifies the path
#  third parameter specifies the file patterns to erase
#  the user is promted before the files are deleted
def QueryDelete(recurse,dir,patn):
    thefiles = ryfiles.listFiles(dir, patn,recurse)
    if len(thefiles)>0:
        for filename in thefiles:
            print(filename)
        if sys.version_info[0] < 3:
            instr=raw_input("Delete these files? (y/n)")
        else:
            instr=input("Delete these files? (y/n)")
           
        if instr=='y':
            for filename in thefiles:
                os.remove(filename)
   

################################################################################
################################################################################
# this is where it all starts
if __name__ == "__main__":

    # process command line agruments
    use = """
    Usage   : python %prog dirName
    
    Example : python %prog UserGuide
    
    The program cleans the specified directory recursively of the following files:
        '*.ps'
        '*.log'
        '*.bbl;comment.cut'
        '*.bbl;*.sav;*.bak;*.synctex;*.log;*.svn'
        '*.blg;*.dfn;*.smb;*.bak;*.aux;*.out;*.lot;*.lof;*.toc;*.tex.bak;*.dvi;*.efc;Backup_of_*.*;*.abr'
     
    If no argument is specified all directories are cleaned.
    
    """

    parser = OptionParser(usage = use)
    options, args = parser.parse_args()

    toclean = os.path.normpath('./') 
    if len(args) > 0:
        toclean = args[0]
    
    print('Cleaning path {0}\n'.format(toclean))

    #we take the conservative approach and do not do blanket erase, 
    #rather do it by type, asking the user first
    QueryDelete(1,toclean, '*.ps;*.fdb_latexmk;*.fls;*synctex.gz')
    QueryDelete(1,toclean, '*.log')
    QueryDelete(1,toclean, '*.bbl;comment.cut')
    QueryDelete(1,toclean, '*.bbl;*.sav;*.bak;*.synctex;*.log;*.svn')
    QueryDelete(1,toclean, '*.blg;*.dfn;*.smb;*.bak;*.aux;*.out;*.lot;*.lof;*.toc;*.tex.bak;*.dvi;*.efc;Backup_of_*.*;*.abr')





