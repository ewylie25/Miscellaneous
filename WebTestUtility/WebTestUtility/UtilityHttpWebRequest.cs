using System;
using System.Net;
using System.Net.Http;
using System.Threading;
using System.Threading.Tasks;

namespace WebTestUtility
{
    public class UtilityHttpWebRequest:IUtilityWebRequest
    {
        private Uri RequestUri { get; }
        private NetworkCredential UserCredentials { get; }

        public UtilityHttpWebRequest(Uri requestUri, NetworkCredential userCredential)
        {
            Log.Information("HTTP Request - creating request");
            this.RequestUri = requestUri;
            this.UserCredentials = userCredential;
        }

        public async Task<string> ExecuteRequest()
        {
            Log.Information("HTTP Request - executing request");
            try
            {
                Log.Information("HTTP Request - creating HttpRequestMessage");
                HttpRequestMessage request = new HttpRequestMessage(HttpMethod.Get, this.RequestUri);

                Log.Information("HTTP Request - creating Handler/client");
                CancellationToken token = new CancellationToken();
                using (var handler = new WebRequestHandler())
                using (HttpClient httpClient = new HttpClient(handler))
                {
                    Log.Information("HTTP Request - setting user credentials");
                    handler.Credentials = this.UserCredentials;

                    Log.Information("HTTP Request - sending request");
                    return await httpClient.SendAsync(request, token).
                        ContinueWith(response
                            =>
                        {
                            try
                            {
                                Log.Information("HTTP Request - processing response");
                                return this.ProcessResponse(response);
                            }
                            catch (AggregateException aggregateException)
                            {
                                Log.Information("HTTP Request - request failed - emitted from task");
                                Log.Exception(aggregateException);
                                return this.ProcessError(aggregateException.ToString());
                            }
                        }, token);

                }

            }
            catch (Exception ex)
            {
                Log.Information("HTTP Request - request failed");
                Log.Exception(ex);
                return this.ProcessError(ex.ToString());
            }
        }

        private string ProcessError(string exceptionString)
        {
            return $"Unhandled .NET error: {exceptionString}";

        }

        private string ProcessResponse(Task<HttpResponseMessage> response)
        {
            if (response.Result.StatusCode == HttpStatusCode.OK ||
                response.Result.StatusCode == HttpStatusCode.Created)
            {
                Log.Information("HTTP Request - success response");
                return this.ProcessSuccess(response.Result.Content.ReadAsStringAsync().TryResult());
            }
            else
            {
                Log.Information("HTTP Request - error response");
                string result = string.Concat("Unsuccessful Response message: \n",
                    $"HttpStatus {response.Result.StatusCode}, " +
                    $"ReasonPhrase {response.Result.ReasonPhrase} \n" +
                    $"Description {response.Result.Content.ReadAsStringAsync().TryResult()}");
                Log.Information(result);
                return result;

            }
        }

        private string ProcessSuccess(object content)
        {
            string result = (string)content;
            Log.Information(result);
            return result;

        }
    }
}