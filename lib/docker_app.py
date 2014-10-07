import os
import sys
import shutil
import json

UROOT = "/data/output/appresults/"
DROOT = "/data/input/samples/"
APPSESSION = "/data/input/TestAppSession.json"

#A base class which handles download, execution and upload
#By inheriting, the child app can supply the details of the app's execution
#and not worry about the plumbing
class DockerApp :

    def run( self ) :
        try :
            self.project_id = sys.argv[1]
            print "In DockerApp downloading files"
            ######  Download  ######
            downloaded_files = self.getDownloadedFiles()
            
            print "In DockerApp, trying to parse AppSession"
            app_session = self.parseAppSession()
            print "In DockerApp, trying to run app"
            ######  Run  ######
            upload_folders = self.runApp( downloaded_files, app_session )

            print "Upload result folders"
            ######  Upload  ######
            self.setupUpload( upload_folders )

        except Exception as e :
            print "Exception: ", type(e), e
            sys.exit(1)

    def getDownloadedFiles(self) :
        downloaded_files = []
        for path, dirs, files in os.walk( DROOT ) :
            for file in files :
                downloaded_files.append( os.path.join( path, file ) )

        return downloaded_files

    def setupUpload( self, upload_folders ) :
        for upload_folder in upload_folders :
            (head, ar_name) = os.path.split( upload_folder )
            target = os.path.join( UROOT, self.project_id )
            shutil.move( upload_folder, os.path.join( target, ar_name ) )
    def parseAppSession( self ) :
        print "In parseAppSession:"
        fh = open(APPSESSION, "r")
        res = json.load(fh)
        if res:
            return res
        else:
            return None
