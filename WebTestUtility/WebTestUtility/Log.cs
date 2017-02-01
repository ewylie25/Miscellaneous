using System;


namespace WebTestUtility
{

    public static class Log
    {
        private static readonly MessageLog MyLog;

        static Log()
        {
            MyLog = new MessageLog();
        }

        public static void Information(string message)
        {
            Message(MsgType.Information, message);
        }

        public static void Error(string message)
        {
            Message( MsgType.Error, message);
        }

        private static void Message(MsgType msgType, string message)
        {
            
            if (null == MyLog)
            {
                return;
            }


            switch (msgType)
            {
                case MsgType.Information:
                    MyLog.Information(message);
                    return;
                case MsgType.Error:
                    MyLog.Error(message);
                    return;
                default:
                    throw new ArgumentOutOfRangeException(nameof(msgType), msgType, null);
            }
        }

        public static void Exception(Exception ex)
        {
            MyLog?.LogUnhandledError(ex);
        }

        private enum MsgType
        {
            Information,
            Error
        }
    }
}