using System;
using System.IO;

namespace WebTestUtility
{
    public class MessageLog
    {
        private readonly string _mChrLogName = "LOG";

        
        private string GetDebugFilePathAndName(DateTime logTime)
        {
            lock (this._lockLog)
            {
                string retValue = $"{this._mChrLogName}_{logTime.Year:0000}_{logTime.Month:00}_{logTime.Day:00}.txt";
                return retValue;
            }
        }

        
        private readonly object _lockLog = new object();

        
        private string TimeStamp
        {
            get
            {
                DateTime dtmNow = DateTime.Now;
                string strDateTime = $"{dtmNow.Hour:00}:{dtmNow.Minute:00}:{dtmNow.Second:00} ";               
                return strDateTime;
            }
        }

        public MessageLog()
        {
            string dir = AppDomain.CurrentDomain.BaseDirectory;


            // Make sure our path ends with "\\"
            if (!dir.EndsWith(Path.DirectorySeparatorChar.ToString()))
            {
                dir += Path.DirectorySeparatorChar.ToString();
            }

            var temp = dir + this._mChrLogName;

            this._mChrLogName = temp;
        }

        public void Information(string strMessage)
        {
            string timeStamp = this.TimeStamp;
            this.Write(timeStamp +" INFO "+ strMessage);
        }


        
        public void Error(string strMessage)
        {
            string timeStamp = this.TimeStamp;
            this.Write(timeStamp + " ERR " + strMessage);
        }

        public void LogUnhandledError(Exception e)
        {
            // Log the error
            this.Error("----------------------------------------------------------------\n");
            this.Error(
                $"UNHANDLED ERROR:\n  Message:\t  {e.Message}\n  Source :\t  {e.Source}\n  Stack:\n{e.StackTrace}\n");



        }
        private void Write(string strMessage)
        {
            try
            {
                lock (this._lockLog)
                {

                    // Get the time
                    DateTime tmTime = DateTime.Now;
                    var chrFile = this.GetDebugFilePathAndName(tmTime);


                    // Open a stream writer for appending
                    using (StreamWriter sw = new StreamWriter(chrFile, true))
                    {
                        sw.WriteLine(strMessage);
                        sw.Close();
                    }

                    
                }

            }
            catch (IOException e)
            {
                Console.WriteLine("Exception: {0}", e.Message);

            }

        }

    }
}