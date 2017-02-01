using System;
using System.Collections.Generic;
using System.Net;
using System.Net.Http;
using System.Security.Cryptography.X509Certificates;
using System.Threading;
using System.Threading.Tasks;

namespace WebTestUtility
{
    public class UtilityHttpsWebRequest:IUtilityWebRequest
    {
        private Uri RequestUri { get; }
        private NetworkCredential UserCredentials { get; }
        private X509Certificate2 Certificate { get; }

        public UtilityHttpsWebRequest(Uri requestUri, NetworkCredential userCredential, X509Certificate2 sslCertificate)
        {
            Log.Information("HTTPS Request - creating request");
            this.RequestUri = requestUri;
            this.UserCredentials = userCredential;
            this.Certificate = sslCertificate;
        }



        public async Task<string> ExecuteRequest()
        {
            Log.Information("HTTPS Request - executing request");
            try
            {
                Log.Information("HTTPS Request - creating HttpRequestMessage");
                HttpRequestMessage request = new HttpRequestMessage(HttpMethod.Get, this.RequestUri);

                Log.Information("HTTPS Request - adding certificate thumbprint");
                List<string> certs = new List<string> {this.Certificate.Thumbprint};
                request.Headers.Add("Thumbprint", certs);

                Log.Information("HTTPS Request - creating Handler/client");
                CancellationToken token = new CancellationToken();
                using (var handler = new WebRequestHandler())
                using (HttpClient httpClient = new HttpClient(handler))
                {
                    Log.Information("HTTPS Request - setting user credentials");
                    handler.Credentials = this.UserCredentials;

                    Log.Information("HTTPS Request - modifying cert validation callback to bypass default validation");
                    handler.ServerCertificateValidationCallback += (sender, certificate, chain, sslPolicyErrors) =>
                    {
                        // investigate certificate parameter
                        X509Certificate2 x509 = new X509Certificate2(certificate);
                        Log.Information($"Certificate expired on: {x509.NotAfter}");
                        return true; // true to bypass, false otherwise
                    };

                    Log.Information("HTTPS Request - sending request");
                    return await httpClient.SendAsync(request, token).
                        ContinueWith(response
                            =>
                        {
                            try
                            {
                                Log.Information("HTTPS Request - processing response");
                                return this.ProcessResponse(response);
                            }
                            catch (AggregateException aggregateException)
                            {
                                Log.Information("HTTPS Request - request failed - emitted from task");
                                Log.Exception(aggregateException);
                                return this.ProcessError(aggregateException.ToString());
                            }
                        }, token);

                }
                
            }
            catch (Exception ex)
            {
                Log.Information("HTTPS Request - request failed");
                Log.Exception(ex);
                return this.ProcessError(ex.ToString());
            }
        }

        private string ProcessError(string exceptionString)
        {
            return $"Unhandled .NET error - see Log: {exceptionString}";

        }

        private string ProcessResponse(Task<HttpResponseMessage> response)
        {
            if (response.Result.StatusCode == HttpStatusCode.OK ||
                response.Result.StatusCode == HttpStatusCode.Created)
            {
                Log.Information("HTTPS Request - success response");
                return this.ProcessSuccess(response.Result.Content.ReadAsStringAsync().TryResult());
            }
            else
            {
                Log.Information("HTTPS Request - error response");
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