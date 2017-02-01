using System;
using System.Threading.Tasks;

namespace WebTestUtility
{
    public static class Extensions
    {
        public static object TryResult<T>(this Task<T> task)
        {
            if (task != null)
            {
                try
                {
                    return task.Result;
                }
                catch (Exception)
                {
                }
            }
            return task;
        }
    }
}